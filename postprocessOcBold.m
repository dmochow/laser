% 10.18.18 running on individually optimally combined me bold

%% NB: segmentAllSubjects must have been run prior to this script
clear all; close all; clc
!export PATH=$PATH:/Users/jacekdmochowski/abin
PATH = getenv('PATH');
setenv('PATH', [PATH ':/Users/jacekdmochowski/abin']);

% constants
subjStrs={'S02','S03','S04','S05','S06','S07','S09','S10','S11','S12','S13'...
    ,'S14','S15','S16','S17','S18','S19','S20','S21','S22'};
nSubjects=numel(subjStrs);
TR=2.8;
fs=1/TR;
CSF_INDEX=1; % segmentation label for CSF
WM_INDEX=3; % segmentation label for white matter
NOMINAL_ONSET_TIME=214; % frame to which align bolds according to (10 mins)

% OPTIONS
CSHIFTOPTION=2; % 1=laserOrigin, 2=ROI centroid
FWHMSTR='8';
FILTER=1;   % lowpass filter at fl Hz
DERIVATIVE=1; % add the derivative of the realignment parameters to regressors
WHITE_MATTER=1; % whether to regress out the white matter PCs
STANDARDIZE=1; % whether to z-score all bold series at the end
TSHIFT=1; % whether to align everyone's laser onset in time
Kwm=3; % number of WM PCs to remove
fl=0.1; % low-pass frequency cutoff
CENTERBOLDS=1;

nPad=30; % number of samples to pad to avoid filter transient (don't change)
[bf,af]=butter(3,fl/(fs/2)); % compute the filter coefficients here
nX=64; nY=76; nZ=60; % dimensions of BOLD
nTR=645; % number of TRs
dummyMask=true(nX,nY,nZ);

optionStr=['-Filter' num2str(FILTER) '-Derivative' num2str(DERIVATIVE) ...
    '-WhiteMatter' num2str(Kwm) ...
    '-STANDARDIZE' num2str(STANDARDIZE) '-TSHIFT' num2str(TSHIFT)  ...
     '-CSHIFT' num2str(CSHIFTOPTION)];

% get the spatial shifts to align bolds
cshifts=computeCshifts(CSHIFTOPTION);

%
alloBolds=cell(nSubjects,1);
muwBold=zeros(nX,nY,nZ,nTR);
alltShifts=zeros(nSubjects,1);
for s=3:numel(subjStrs)
    s
    subjIndx=subjStrs{s};
    basePath=['../data/' subjIndx '/NII/'];
    inputBoldFilename=fullfile(basePath,'oc_dsbold+tlrc.BRIK');
    brainMaskFilename=fullfile(basePath,'resampled_brain_mask+tlrc.BRIK');
    motionFilename=fullfile(basePath,'dsbold_e2_vr_motion.1D');
    resampledSegFilename=fullfile(basePath,'resampled_Classes+tlrc.BRIK');
    biopacPath=['../data/' subjIndx '/BIOPAC'];
    biopacFilename='biopac.mat';
    
    % biopac
    load(fullfile(biopacPath,biopacFilename),'data');
    [~,maxindsample]=max(diff(data(:,1)));
    onsetTimeSec=maxindsample/1000;
    onsetTimeTR=round(onsetTimeSec/TR);
    tshift=NOMINAL_ONSET_TIME-onsetTimeTR;
    alltShifts(s)=tshift;
    
    [~, inputBold, iInfo, ~] = BrikLoad (inputBoldFilename);
    iwBoldTs=vol2ts(inputBold);
    
    [~, brainMask, ~, ~] = BrikLoad (brainMaskFilename);
    brainMask=logical(brainMask);
    
    % create the motion regressors
    fid=fopen(motionFilename);
    mals=textscan(fid,'%f %f %f %f %f %f');
    mals=cell2mat(mals);
    if DERIVATIVE
        dmals=cat(1,zeros(1,6),diff(mals));
        mals=cat(2,mals,dmals);
    end
    mals=mals-repmat(mean(mals),size(mals,1),1);
    
    % regress out motion
    owBoldTs = regressOut(iwBoldTs.',mals.',1).';
    % put back to 4D
    outputBold=ts2vol(owBoldTs,dummyMask);
    
    if WHITE_MATTER
        [~,seg,~,~] = BrikLoad (resampledSegFilename);
        wmMask=seg==WM_INDEX;
        wmBoldTs=vol2ts(outputBold,wmMask);
        [U,S,V]=svd(wmBoldTs,0);
        brainBoldTs=vol2ts(outputBold,brainMask); % only regress in the brain mask
        obrainBoldTs = regressOut(brainBoldTs.',U(:,1:Kwm).',1).';
    end
      
    if FILTER
        tmp=cat(1,repmat(obrainBoldTs(1,:),nPad,1),obrainBoldTs);
        tmpOut=filter(bf,af,tmp,[],1);
        obrainBoldTs=tmpOut(nPad+1:end,:);
    end
    
    if STANDARDIZE
        obrainBoldTs=zscore(obrainBoldTs);
    end
    
    if TSHIFT
        % temporally shift bold to align everyone's onset time
        obrainBoldTs=circshift(obrainBoldTs,[tshift 0]);
    end
    final4D=ts2vol(obrainBoldTs,brainMask);
    
    % write out the post-processed bold here
    [err,ErrMessage,Info]=WriteBrikWrap(basePath,final4D,iInfo,['ocBold' optionStr],'tlrc');
%     
    % here we center the bold so that everyone's ROI lines up
    % this does not take into account different laser orientations
    if CENTERBOLDS
        final4Ds=circshift(final4D,[cshifts(1,s) cshifts(2,s) cshifts(3,s) 0]);
        muwBold=muwBold+final4Ds/nSubjects;
    else
        muwBold=muwBold+final4D/nSubjects;
    end

end


%%
% write out the average bold here
[err,ErrMessage,info]=WriteBrikWrap('../data/SAVG/NII/',muwBold,iInfo,['muocBold' optionStr],'tlrc');

%%
% smooth the average bold here
for f=1:numel(FWHMSTR)
    try
    tFWHMSTR=FWHMSTR{f};
    catch
        tFWHMSTR=FWHMSTR; % just one, not a cell
    end
    origPath=pwd;
    cd('../data/SAVG/NII/');
    muwPrefix=['s' tFWHMSTR 'muocBold' optionStr];
    str=['!3dBlurInMask -prefix ' muwPrefix ' -input muocBold' optionStr '+tlrc -FWHM ' tFWHMSTR ' -mask muocBold' optionStr '+tlrc'];
    eval(str);
    cd(origPath);
end
