% 04.19.19
% take the anatomy aligned and talairached bold
% apply minimal preprocessing: regressing out motion regressors
% convert to signal change and save

% 04.25.19
% something

% 04.26.19
% trying to replicate result of processSubject.m

% 04.27.19
% AFNI does this:
% despike tshift align tlrc volreg blur mask        
% scale regress

% 04.29.19
% making this the general file to run after processSubject block 1
% input to this file is: dsbold_e*_tlrc_al+tlrc.BRIK
% this needs to handle both one-shot and toggled data
clear all; close all; clc
addpath(genpath('/Users/jacekdmochowski/PROJECTS/COMMON'));
!export PATH=$PATH:/Users/jacekdmochowski/abin
PATH = getenv('PATH');
setenv('PATH', [PATH ':/Users/jacekdmochowski/abin']);
origPath=pwd;

%% PREAMBLE
subjStrs={'S02','S03','S04','S05','S06','S07','S09','S10','S11','S12','S13'...
    ,'S14','S15','S16','S17','S18','S19','S20','S21','S22'};
%subjStrs={'S23','S24','S25','S26','S27','S28','S29','S30'};
nSubjects=numel(subjStrs);
TR=2.8;
fs=1/TR;
nX=64; nY=76; nZ=60; % dimensions of BOLD
nTR=645; % number of TRs
dummyMask=true(nX,nY,nZ); % used for going from 2D to 4D
radMask=21.3; depthMask=39.2; % 99 percent min area
bslInds=1:214; % for computing percent change
gmIndex=2; % segmentation label for grey matter
wmIndex=3; % segmentation label for white matter
Kwm=3; % number of WM PCs to regress out
fl=0.1; % low-pass frequency cutoff
nPad=30; % number of samples to pad to avoid filter transient (don't change)
[bf,af]=butter(3,fl/(fs/2)); % compute the filter coefficients here
fwhm=8;

%% main loop
for s=1:nSubjects
    subjStr=subjStrs{s};
    pathToData=['../data/' subjStr '/NII/'];
    %pathToData=['../data/' subjStr];
    
    % load in brain mask
    brainMaskFilename=fullfile(pathToData,'resampled_brain_mask+tlrc.BRIK');
    % note that we are working with functional data in the talairach space here
    [~, brainMask, ~, ~] = BrikLoad (brainMaskFilename);
    brainMask=logical(brainMask);
    
    % load in segmentation results
    resampledSegFilename=fullfile(pathToData,'resampled_Classes+tlrc.BRIK');
    [~,seg,~,~] = BrikLoad(resampledSegFilename);
    wmMask=seg==wmIndex;
    gmMask=seg==gmIndex;
    
    %%
    % load input bold ts
    for e=1:3
        [s,e]

        inputBoldFilename=fullfile(pathToData,['dsbold_e' num2str(e) '_tlrc_al+tlrc.BRIK']);
        [~, inputBold, iInfo, ~] = BrikLoad (inputBoldFilename);
        brainTs=vol2ts(inputBold,brainMask);
        
        % create the motion regressors
        motionFilename=fullfile(pathToData,['dsbold_e' num2str(e) '_vr_motion.1D']);
        fid=fopen(motionFilename);
        mals=textscan(fid,'%f %f %f %f %f %f');
        mals=cell2mat(mals);
        dmals=cat(1,zeros(1,6),diff(mals));
        mals=cat(2,mals,dmals);
        mals=mals-repmat(mean(mals),size(mals,1),1);
        
        % regress out motion
        obrainTs = regressOut(brainTs.',mals.',1).';
        
        % construct WM regressor
        tmp=ts2vol(obrainTs,brainMask); % back to 4D
        wmBoldTs=vol2ts(tmp,wmMask);
        [U,S,V]=svd(wmBoldTs,0);
        
        % regress out white matter signal 
        oobrainTs = regressOut(obrainTs.',U(:,1:Kwm).',1).';
        
        % filter
        tmp=cat(1,repmat(oobrainTs(1,:),nPad,1),oobrainTs);
        tmpOut=filter(bf,af,tmp,[],1);
        foobrainTs=tmpOut(nPad+1:end,:);
        
        % z-score
        zfoobrainTs=zscore(foobrainTs);
        
        % 2D --> 4D
        final4D=ts2vol(zfoobrainTs,brainMask);
    
        % write out the post-processed bold here
        [err,ErrMessage,Info]=WriteBrikWrap(pathToData,final4D,iInfo,['unsmoothedOutputBoldEcho' num2str(e)],'tlrc');
    
        % blur
        inputStr=['unsmoothedOutputBoldEcho' num2str(e) '+tlrc '];
        smPrefix=['smoothedOutputBoldEcho' num2str(e)];
        cd(pathToData);
        stro=['delete smoothedOutputBoldEcho' num2str(e) '+tlrc.BRIK']; eval(stro);
        stro=['delete smoothedOutputBoldEcho' num2str(e) '+tlrc.HEAD']; eval(stro);
        str=['!3dBlurInMask -prefix ' smPrefix ' -input ' inputStr '-FWHM ' num2str(fwhm) ' -mask ' inputStr];
        eval(str);
        cd(origPath);
    end
    
    
end

