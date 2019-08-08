%% adapt:
% (1) handle all subjects 2-32
% (2) put output data into separate folder

%% NB: this code was written for Study II: toggling of the laser
%% it is not applicable to Study I in the present form
clear all; close all; clc
addpath(genpath('/Users/jacekdmochowski/PROJECTS/COMMON'));
!export PATH=$PATH:/Users/jacekdmochowski/abin
PATH = getenv('PATH');
setenv('PATH', [PATH ':/Users/jacekdmochowski/abin']);

%% PREAMBLE
subjStr='S24';
pathToData=['../data/' subjStr];
origPath=pwd;
laserOriginPrefix='laserOriginV';
anatFilename=fullfile(pathToData,'anat+orig');
anatWithSkullFilename=fullfile(pathToData,'anatWithSkull+orig');
brainMaskFilename=fullfile(pathToData,'resampled_brain_mask+tlrc.BRIK'); % note that we are working with functional data in the talairach space here
resampledSegFilename=fullfile(pathToData,'resampled_Classes+tlrc.BRIK');
anatBrainMaskFilename=fullfile(pathToData,'brain_mask+orig'); 
biopacFilename=fullfile(pathToData,'subj0024_032819b.mat');

TR=2.8;
fs=1/TR;
nX=64; nY=76; nZ=60; % dimensions of BOLD
nTR=645; % number of TRs
dummyMask=true(nX,nY,nZ); % used for going from 2D to 4D
Kwm=3; % number of WM PCs to remove
fl=0.1; % low-pass frequency cutoff
nPad=30; % number of samples to pad to avoid filter transient (don't change)
[bf,af]=butter(3,fl/(fs/2)); % compute the filter coefficients here
WM_INDEX=3; % segmentation label for white matter
radMask=21.3; depthMask=39.2; % 99 percent min area

%% create output directory, copy files, and begin
outputDir=['../data/myoutput/' subjStrs{s}]

%% Part I: brain mask generation
cd(pathToData);
if ~exist('brain_mask+orig.HEAD','file')
    str1='!3dcopy anat.nii anat';
    str1b='!3dcopy anat.nii anatWithSkull';
    str2='!3dSkullStrip -input anat+orig';
    str3='!3dAutomask -prefix brain_mask skull_strip_out+orig';
    str3b='!gunzip brain_mask+orig.BRIK.gz';
    str3c='!3dSeg -anat skull_strip_out+orig -mask AUTO -classes ''CSF ; GM ; WM'' -bias_classes ''GM ; WM'' -bias_fwhm 25 -mixfrac UNI -main_N 5 -blur_meth BFT';
    eval(str1);
    eval(str1b);
    eval(str2);
    eval(str3);
    eval(str3b);
    eval(str3c);
end


%% Part II: bold-anatomy registration, atlasing, and resampling
str4='!3dcopy bold_e1.nii bold_e1'; eval(str4);
str5='!3dcopy bold_e2.nii bold_e2'; eval(str5);
str6='!3dcopy bold_e3.nii bold_e3'; eval(str6);
str7='!3ddespike -NEW -localedit -prefix dsbold_e1 bold_e1+orig'; eval(str7);
str8='!3ddespike -NEW -localedit -prefix dsbold_e2 bold_e2+orig'; eval(str8);
str9='!3ddespike -NEW -localedit -prefix dsbold_e3 bold_e3+orig'; eval(str9);
% align the anatomical to the talairach atlas
str10='!@auto_tlrc -base TT_N27+tlrc -input anat+orig'; eval(str10);
% align the brain mask to the talairach atlas
str10b='!@auto_tlrc -apar anat+tlrc -input brain_mask+orig'; eval(str10b);
str10c='!cp ./Segsy/Classes+orig.BRIK.gz .'; eval(str10c);  
str10cc='!cp ./Segsy/Classes+orig.HEAD .'; eval(str10cc);  
str10d='!gunzip Classes+orig.BRIK.gz'; eval(str10d); 
str10e='!@auto_tlrc -apar anat+tlrc -input Classes+orig.BRIK'; eval(str10e);
% now align the epis to the anatomy
str12='!align_epi_anat.py -anat anat+orig -epi dsbold_e2+orig -epi_base 5 -child_epi dsbold_e1+orig dsbold_e3+orig -epi2anat -suffix _al -tlrc_apar anat+tlrc';
eval(str12);
% resample the brain mask on the bold grid
str14b='!3dresample -master dsbold_e2_tlrc_al+tlrc -prefix resampled_brain_mask -input brain_mask+tlrc'; eval(str14b);
str14c='!3dresample -master dsbold_e2_tlrc_al+tlrc -prefix resampled_Classes -input Classes+tlrc'; eval(str14c);
cd(origPath);


%% DENOISING
% MATLAB begins here
for e=1:3
    
    % load input bold ts
    inputBoldFilename=fullfile(pathToData,['dsbold_e' num2str(e) '_tlrc_al+tlrc.BRIK']);
    [~, inputBold, iInfo, ~] = BrikLoad (inputBoldFilename);
    iwBoldTs=vol2ts(inputBold);
    
    % load in brain mask
    [~, brainMask, ~, ~] = BrikLoad (brainMaskFilename);
    brainMask=logical(brainMask);
    
    % create the motion regressors
    motionFilename=fullfile(pathToData,['dsbold_e' num2str(e) '_vr_motion.1D']);
    fid=fopen(motionFilename);
    mals=textscan(fid,'%f %f %f %f %f %f');
    mals=cell2mat(mals);
    dmals=cat(1,zeros(1,6),diff(mals));
    mals=cat(2,mals,dmals);
    mals=mals-repmat(mean(mals),size(mals,1),1);
    
    % regress out motion
    owBoldTs = regressOut(iwBoldTs.',mals.',1).';
    outputBold=ts2vol(owBoldTs,dummyMask);

    % regress out white matter
    [~,seg,~,~] = BrikLoad(resampledSegFilename);
    wmMask=seg==WM_INDEX;
    wmBoldTs=vol2ts(outputBold,wmMask);
    [U,S,V]=svd(wmBoldTs,0);
    brainBoldTs=vol2ts(outputBold,brainMask); % only regress in the brain mask
    obrainBoldTs = regressOut(brainBoldTs.',U(:,1:Kwm).',1).';
    
    % filter
    tmp=cat(1,repmat(obrainBoldTs(1,:),nPad,1),obrainBoldTs);
    tmpOut=filter(bf,af,tmp,[],1);
    obrainBoldTs=tmpOut(nPad+1:end,:);
    
    % z-score
    obrainBoldTs=zscore(obrainBoldTs);
    
    % 2D --> 4D
    final4D=ts2vol(obrainBoldTs,brainMask);
    
    % write out the post-processed bold here
    [err,ErrMessage,Info]=WriteBrikWrap(pathToData,final4D,iInfo,['smoothedOutputBoldEcho' num2str(e)],'tlrc');
    
    % smooth
    inputStr=['smoothedOutputBoldEcho' num2str(e) '+tlrc '];
    smPrefix=['smoothed8mmOutputBoldEcho' num2str(e)];
    cd(pathToData);
    str=['!3dBlurInMask -prefix ' smPrefix ' -input ' inputStr '-FWHM 8 -mask ' inputStr];
    eval(str);
    cd(origPath);
 
end

%% ROI GENERATION
[urf,ulf,urb,lrf,midsag]=getMarkerCoords(subjStr);
laserOrigin = getLaserOrigin(ulf,urf,urb);
laserOrigin = fineTuneOrigin(anatWithSkullFilename,laserOrigin,ulf,urf,lrf);
save(fullfile(pathToData,'laserOrigin.mat'),'laserOrigin');
roiMask = projectLaser(anatBrainMaskFilename,laserOrigin,urf,ulf,urb,lrf,radMask,depthMask);  % project the real laser into the mri
[~, ~, Info, ~] = BrikLoad (anatFilename); % only to get the 'Info' field
roiPrefix=['roi_r' num2str(radMask,'%0.0f') '_z' num2str(depthMask,'%0.0f')  ];
WriteBrikWrap(pathToData,roiMask,Info,roiPrefix,'orig');

cd(pathToData);
str11=['!@auto_tlrc -apar anat+tlrc -input ' roiPrefix '+orig'] ; eval(str11); % register roi to anatomy
str11a=['!gunzip ' roiPrefix '+tlrc.BRIK.gz']; eval(str11a); % weird that this is needed
str12=['!3dresample -master dsbold_e2_tlrc_al+tlrc -prefix resampled_' roiPrefix ' -input ' roiPrefix '+tlrc']; eval(str12);  % resample roi on bold grid
cd(origPath);

%% define laser time course regressor from BIOPAC
load(biopacFilename,'data');
thresh=0.75; % anything over this voltage is classified "laser on"
laserTc=data;
biopacTime=(0:size(data,1)-1)*0.001; % biopac time in ms but check by loading in all the variables
boldTime=(0:644)*TR; % note that the biopac time is likely shorter than bold
laserTcDown=interp1(biopacTime,laserTc(:,1),boldTime);
regressor=laserTcDown>thresh;
save(fullfile(pathToData,'regressor.mat'),'regressor');