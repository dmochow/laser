%% 07.19.19
% preprocess ASRC MR Thermo output after converting from dicom to nifti
%% 07.22.19
% removing Talairaching (for now)
clear all; close all; clc
addpath(genpath('/Users/jacekdmochowski/PROJECTS/COMMON'));
!export PATH=$PATH:/Users/jacekdmochowski/abin
PATH = getenv('PATH');
setenv('PATH', [PATH ':/Users/jacekdmochowski/abin']);

%% PREAMBLE
subjStrs={'tS01a'};
radMask=21.3; depthMask=39.2; % 99 percent min area
roiPrefix=['roi_r' num2str(radMask,'%0.0f') '_z' num2str(depthMask,'%0.0f')  ];
origPath=pwd;


%% main loop
nSubjects=numel(subjStrs);
for s=1:nSubjects
    subjStr=subjStrs{s};
    
    % create output directory, copy files, and begin
    pathToScannerData=['/Users/jacekdmochowski/PROJECTS/LLLT/data/thermo/scanner_data/' subjStr];
    outputDir=['/Users/jacekdmochowski/PROJECTS/LLLT/data/thermo/output/' subjStr];
    mkdir(outputDir);
    cd(outputDir);
    copyfile(fullfile(pathToScannerData,'anat.nii'),fullfile(outputDir,'anat.nii'));
    copyfile(fullfile(pathToScannerData,'mag.nii'),fullfile(outputDir,'mag.nii'));
    copyfile(fullfile(pathToScannerData,'ph.nii'),fullfile(outputDir,'ph.nii'));
    
    % define necessary data structures here
    anatFilename=fullfile(outputDir,'anat+orig');
    anatWithSkullFilename=fullfile(outputDir,'anatWithSkull+orig');
    brainMaskFilename=fullfile(outputDir,'resampled_brain_mask+tlrc.BRIK'); % note that we are working with functional data in the talairach space here
    resampledSegFilename=fullfile(outputDir,'resampled_Classes+tlrc.BRIK');
    anatBrainMaskFilename=fullfile(outputDir,'brain_mask+orig');
    
    
    % Part I: brain mask generation
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
    
    cd(origPath);
    % Part II: roi generation
    [urf,ulf,urb,lrf,midsag]=getMarkerCoords(subjStr);
    if ~isempty(urf)
        laserOrigin = getLaserOrigin(ulf,urf,urb);
        %try
        laserOrigin = fineTuneOrigin(anatWithSkullFilename,laserOrigin,ulf,urf,lrf);
        %catch
        %end
        save(fullfile(pathToScannerData,'laserOrigin.mat'),'laserOrigin');
        roiMask = projectLaser(anatBrainMaskFilename,laserOrigin,urf,ulf,urb,lrf,radMask,depthMask);  % project the real laser into the mri
        [~, ~, Info, ~] = BrikLoad (anatFilename); % only to get the 'Info' field
        WriteBrikWrap(outputDir,roiMask,Info,roiPrefix,'orig');
    else
        error('Vitamin E markers missing');
    end
    
    cd(outputDir);
    % Part III: bold-anatomy registration, atlasing, and resampling
    str4='!3dcopy mag.nii mag'; eval(str4);
    str5='!3dcopy ph.nii ph'; eval(str5);
    
    % align the anatomical to the talairach atlas
    %str10='!@auto_tlrc -base TT_N27+tlrc -input anat+orig'; eval(str10);
    % align the brain mask to the talairach atlas
    %str10b='!@auto_tlrc -apar anat+tlrc -input brain_mask+orig'; eval(str10b);
    str10c='!cp ./Segsy/Classes+orig.BRIK.gz .'; eval(str10c);
    str10cc='!cp ./Segsy/Classes+orig.HEAD .'; eval(str10cc);
    str10d='!gunzip Classes+orig.BRIK.gz'; eval(str10d);
    %str10e='!@auto_tlrc -apar anat+tlrc -input Classes+orig.BRIK'; eval(str10e);
    % now align the epis to the anatomy
    str12=['!align_epi_anat.py -anat anat+orig -epi mag+orig -epi_base mean -child_epi ph+orig -epi2anat -suffix _al'];
    eval(str12);
    % resample the brain mask on the bold grid
    str14b='!3dresample -master mag_al+orig -prefix resampled_brain_mask -input brain_mask+orig'; eval(str14b);
    str14c='!3dresample -master mag_al+orig -prefix resampled_Classes -input Classes+orig'; eval(str14c);
    str14d=['!3dresample -master mag_al+orig -prefix resampled_' roiPrefix ' -input ' roiPrefix '+orig']; eval(str14d);  % resample roi on bold grid eval(str14d);
    cd(origPath);
    
    cd(outputDir);
    
    
    % load input bold ts
    inputBoldFilename=fullfile(outputDir,['dsbold_e' num2str(e) '_tlrc_al+tlrc.BRIK']);
    [~, inputBold, iInfo, ~] = BrikLoad (inputBoldFilename);
    iwBoldTs=vol2ts(inputBold);
    
    % load in brain mask
    [~, brainMask, ~, ~] = BrikLoad (brainMaskFilename);
    brainMask=logical(brainMask);
    
    % create the motion regressors
    motionFilename=fullfile(outputDir,['dsbold_e' num2str(e) '_vr_motion.1D']);
    fid=fopen(motionFilename);
    mals=textscan(fid,'%f %f %f %f %f %f');
    mals=cell2mat(mals);
    dmals=cat(1,zeros(1,6),diff(mals));
    mals=cat(2,mals,dmals);
    mals=mals-repmat(mean(mals),size(mals,1),1);
    
    % regress out motion
    owBoldTs = regressOut(iwBoldTs.',mals.',1).';
    outputBold=ts2vol(owBoldTs,dummyMask);
    
    %     str11=['!@auto_tlrc -apar anat+tlrc -input ' roiPrefix '+orig'] ; eval(str11); % register roi to anatomy
    %     str11a=['!gunzip ' roiPrefix '+tlrc.BRIK.gz']; eval(str11a); % weird that this is needed
    %     str12=['!3dresample -master dsmag_tlrc_al+tlrc -prefix resampled_' roiPrefix ' -input ' roiPrefix '+tlrc']; eval(str12);  % resample roi on bold grid
    %
    %     % new July 2019
    %     str13=['!3dresample -master dsmag_al+orig -prefix resampled_' roiPrefix ' -input ' roiPrefix '+orig']; eval(str13);  % resample roi on bold grid
    %     str14=['!3dresample -master dsmag_al+orig -prefix resampled_Classes -input Classes+orig']; eval(str14);  % resample roi on bold grid
    %
    
    
    %cd(outputDir); delete bold_e1.nii; delete bold_e2.nii; delete bold_e3.nii; delete anat.nii;
    
    cd(origPath);
end

