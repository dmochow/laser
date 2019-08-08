% 05/11/18
% this script picks up where AfNI (called in preprocessBold.m) left off,
% and does the rest of the analysis in matlab
%
% writing this to test out different processing chains and their effects on
% motion artifact reduction
clear all; close all; clc
!export PATH=$PATH:/Users/jacekdmochowski/abin
PATH = getenv('PATH');
setenv('PATH', [PATH ':/Users/jacekdmochowski/abin']);

% constants
subjStrs={'S02','S03','S04','S05','S06','S07','S09','S10','S11','S12','S13','S14','S15','S16','S17','S18','S19'};
nSubjects=numel(subjStrs);
TR=2.8;
fs=1/TR;
WM_INDEX=3; % segmentation label for white matter
CSF_INDEX=1; % segmentation label for CSF

% OPTIONS
FILTER=1;   % lowpass filter at fl Hz
DERIVATIVE=1; % add the derivative of the realignment parameters to regressors
WHITE_MATTER=1; % whether to regress out the white matter PCs
CSF=0; % whether to regress out the CSF PCs
Kwm=5; % number of WM PCs to remove
Kcsf=5; % number of CSF PCs to remove
fl=0.1; % low-pass frequency cutoff
nPad=30; % number of samples to pad to avoid filter transient (don't change)
[bf,af]=butter(3,fl/(fs/2)); % compute the filter coefficients here
NX=64; NY=76; NZ=60;
dummyMask=logical(ones(NX,NY,NZ));

outFilenameStr=['../data/precomputed/allBoldsRoi-Filter' num2str(FILTER) '-Derivative' num2str(DERIVATIVE) ...
    '-WhiteMatter' num2str(WHITE_MATTER) '-CSF' num2str(CSF) '-' date];

alloBolds=cell(nSubjects,1);
muwBold=zeros(NX,NY,NZ);
for s=1:numel(subjStrs)
    s
    subjIndx=subjStrs{s};
    basePath=['../data/' subjIndx '/NII/'];
    inputBoldFilename=fullfile(basePath,'dsbold_e1_tlrc_al+tlrc.BRIK');
    brainMaskFilename=fullfile(basePath,'resampled_brain_mask+tlrc.BRIK');
    muRoiMaskFilename=fullfile(basePath,'muroi+tlrc.BRIK');
    roiMaskFilename=fullfile(basePath,'resampled_roi+tlrc.BRIK');
    motionFilename=fullfile(basePath,'dsbold_e1_vr_motion.1D');
    resampledSegFilename=fullfile(basePath,'resampled_Classes+tlrc.BRIK');
    
    [~, inputBold, ~, ~] = BrikLoad (inputBoldFilename);
    [~,~,~,nTR]=size(inputBold);
    iwBoldTs=vol2ts(inputBold);
    
    [~, roiMask, ~, ~] = BrikLoad (roiMaskFilename);
    roiMask=logical(roiMask);
    
    if isempty(roiMask)
        [~, roiMask, ~, ~] = BrikLoad (muRoiMaskFilename);
        roiMask=logical(roiMask);
    end
    
    [~, brainMask, ~, ~] = BrikLoad (brainMaskFilename);
    brainMask=logical(brainMask);
    
    finalMask=roiMask>0 & brainMask>0;

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
        % processing myOutputBold
        [~,seg,~,~] = BrikLoad (resampledSegFilename);
        wmMask=seg==WM_INDEX;
        wmBoldTs=vol2ts(outputBold,wmMask);
        [U,S,V]=svd(wmBoldTs,0);
        owBoldTs = regressOut(owBoldTs.',U(:,1:Kwm).',1).';
        % regress out first K PCs
    end
    
   
    if CSF
        outputBold=ts2vol(owBoldTs,dummyMask); % put back to full 4D
        csfMask=seg==CSF_INDEX;
        csfBoldTs=vol2ts(outputBold,csfMask);
        [U,S,V]=svd(csfBoldTs,0);
        owBoldTs = regressOut(owBoldTs.',U(:,1:Kcsf).',1).';
    end

    if FILTER
        %tmp=cat(1,repmat(oBoldTs(1,:),nPad,1),oBoldTs);
        %tmpOut=filter(bf,af,tmp,[],1);
        %foBoldTs=tmpOut(nPad+1:end,:);
        tmp=cat(1,repmat(owBoldTs(1,:),nPad,1),owBoldTs);
        tmpOut=filter(bf,af,tmp,[],1);
        owBoldTs=tmpOut(nPad+1:end,:);
    end
    
    final4D=ts2vol(owBoldTs,dummyMask);
    
    %
    
    muwBold=muwBold+final4D/nSubjects;
    
    oBoldTs=vol2ts(final4D,finalMask);
    alloBolds{s}=oBoldTs;
end

%%
allBolds2D=cat(2,alloBolds{:});
figure;
subplot(211);
time=(0:nTR-1)*TR;
plot(time,mean(allBolds2D,2),'k');

print('-dpng',['../figures/tmp-' date]);

save(outFilenameStr,'alloBolds','fl','nPad','bf','af','Kwm','Kcsf','FILTER','DERIVATIVE','WHITE_MATTER','CSF','TR','fs');

%%
% write out the average bold here
[~, inputBold, info, ~] = BrikLoad (inputBoldFilename); % just to get "Info"
[err,ErrMessage,info]=WriteBrikWrap('../data/SAVG/NII/',muwBold,info,'muwBold','tlrc');