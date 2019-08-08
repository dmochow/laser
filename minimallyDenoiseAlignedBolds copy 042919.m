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
clear all; close all; clc
addpath(genpath('/Users/jacekdmochowski/PROJECTS/COMMON'));
!export PATH=$PATH:/Users/jacekdmochowski/abin
PATH = getenv('PATH');
setenv('PATH', [PATH ':/Users/jacekdmochowski/abin']);
origPath=pwd;

%% PREAMBLE
%subjStrs={'S02','S03','S04','S05','S06','S07','S09','S10','S11','S12','S13'...
%    ,'S14','S15','S16','S17','S18','S19','S20','S21','S22'};
subjStrs={'S23','S24','S25','S26','S27','S28','S29','S30'};
nSubjects=numel(subjStrs);
TR=2.8;
fs=1/TR;
nX=64; nY=76; nZ=60; % dimensions of BOLD
nTR=645; % number of TRs
dummyMask=true(nX,nY,nZ); % used for going from 2D to 4D
radMask=21.3; depthMask=39.2; % 99 percent min area
bslInds=1:214; % for computing percent change
wmIndex=3; % segmentation label for white matter
Kwm=3; % number of WM PCs to regress out
fl=0.1; % low-pass frequency cutoff
nPad=30; % number of samples to pad to avoid filter transient (don't change)
[bf,af]=butter(3,fl/(fs/2)); % compute the filter coefficients here
fwhm=4;

%% main loop
for s=1:nSubjects
    subjStr=subjStrs{s};
    %pathToData=['../data/' subjStr '/NII/'];
    pathToData=['../data/' subjStr];
    
    % load in brain mask
    brainMaskFilename=fullfile(pathToData,'resampled_brain_mask+tlrc.BRIK');
    % note that we are working with functional data in the talairach space here
    [~, brainMask, ~, ~] = BrikLoad (brainMaskFilename);
    brainMask=logical(brainMask);
    
    % load in segmentation results
    resampledSegFilename=fullfile(pathToData,'resampled_Classes+tlrc.BRIK');
    [~,seg,~,~] = BrikLoad(resampledSegFilename);
    wmMask=seg==wmIndex;
    
    %%
    % load input bold ts
    for e=1:3
        [s,e]
             
        % blur
        smBold = smoothBold(pathToData,['dsbold_e' num2str(e) '_tlrc_al+tlrc.BRIK'],fwhm,0);
        brainTs=vol2ts(smBold,brainMask);
        
        % afni style scaling
        mus=mean(brainTs,1);
        brainTs=brainTs./repmat(mus,size(brainTs,1),1);
        brainTs=brainTs*100; % afni style
        
        inputBoldFilename=fullfile(pathToData,['dsbold_e' num2str(e) '_tlrc_al+tlrc.BRIK']);
        [~, ~, iInfo, ~] = BrikLoad (inputBoldFilename);
        
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
        wmBoldTs(isnan(wmBoldTs))=0;
        [U,S,V]=svd(wmBoldTs,0);
        
        % regress out white matter signal 
        obrainTs = regressOut(obrainTs.',U(:,1:Kwm).',1).';
        
        % filter
        tmp=cat(1,repmat(obrainTs(1,:),nPad,1),obrainTs);
        tmpOut=filter(bf,af,tmp,[],1);
        obrainTs=tmpOut(nPad+1:end,:);
        
        % 2D --> 4D
        obrain4D=ts2vol(obrainTs,brainMask);
       
        % write out the unsmoothed bold here
        [err,ErrMessage,Info]=WriteBrikWrap(pathToData,obrain4D,iInfo,['smoothedOutputBoldEcho' num2str(e)],'tlrc');
        
%         % blur
%         obrain4D = smoothBold(pathToData,['unsmoothedOutputBoldEcho' num2str(e) '+tlrc'],fwhm,0);
%         
%         % 4D --> 2D
%         obraints=vol2ts(obrain4D,brainMask);
%         
% 
%         
%         % z-score
%         %nfoobrainTs=zscore(foobrainTs);
%         
%         % afni style scaling
%         mus=mean(obrainTs,1);
%         obrainTs=obrainTs./repmat(mus,size(obrainTs,1),1);
%         obrainTs=obrainTs*100; % afni style
%         
%         % 2D --> 4D
%         final4D=ts2vol(obrainTs,brainMask);
%         
%         [err,ErrMessage,Info]=WriteBrikWrap(pathToData,final4D,iInfo,['smoothedOutputBoldEcho' num2str(e)],'tlrc');
                
    end
    
    
end

