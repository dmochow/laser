% 08/28/18: may have to rerun but aligning with laserOrigin or absorption?
% 06/05/18: WM regression now only done in brain mask as doing it in the
% whole volume was leading to massive blurring of the bold
%
% 06/02/18: correct biopac files inserted in all subfolders
%
% 05/31/18:
% all bold time series have been standardized
% echo number and FWHM are now parameters
%
% 05/30/18: check the effect of potential different baselines on averaging process
%
% 05/11/18
% this script picks up where AfNI (called in preprocessBold.m) left off,
% and does the rest of the analysis in matlab
%
% writing this to test out different processing chains and their effects on
% motion artifact reduction

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
ECHOSTR='1';
FWHMSTR={'8','10','20'};
FILTER=1;   % lowpass filter at fl Hz
DERIVATIVE=1; % add the derivative of the realignment parameters to regressors
WHITE_MATTER=1; % whether to regress out the white matter PCs
CSF=0; % whether to regress out the CSF PCs
STANDARDIZE=1; % whether to z-score all bold series at the end
TSHIFT=1; % whether to align everyone's laser onset in time
Kwm=3; % number of WM PCs to remove
Kcsf=1; % number of CSF PCs to remove
fl=0.1; % low-pass frequency cutoff
CENTERBOLDS=1;

% ROI options
%radMask=19.2; depthMask=26.5; % original values
radMask=21.3; depthMask=39.2; % original values
roiPrefix=['roi_r' num2str(radMask,'%0.0f') '_z' num2str(depthMask,'%0.0f')  ];
controlRoiPrefix=['control_roi_r' num2str(radMask,'%0.0f') '_z' num2str(depthMask,'%0.0f')  ];
roiFilename=[roiPrefix '+tlrc'];
controlRoiFilename=[controlRoiPrefix '+tlrc'];
muRoiPrefix=['mu_' roiPrefix];
muControlRoiPrefix=['mu_' controlRoiPrefix];

nPad=30; % number of samples to pad to avoid filter transient (don't change)
[bf,af]=butter(3,fl/(fs/2)); % compute the filter coefficients here
nX=64; nY=76; nZ=60; % dimensions of BOLD
nTR=645; % number of TRs
dummyMask=true(nX,nY,nZ);

outFilenameStr=['../data/precomputed/allBoldsRoi-Filter' num2str(FILTER) '-' roiPrefix '-Derivative' num2str(DERIVATIVE) ...
    '-WhiteMatter' num2str(WHITE_MATTER) '-KWM' num2str(Kwm) '-CSF' num2str(CSF) '-KCSF' num2str(Kcsf)  '-STANDARDIZE' num2str(STANDARDIZE) '-TSHIFT' num2str(TSHIFT) '-Echo' ECHOSTR '-' date]

figFilenameStr=['../figures/grandMean-Filter' num2str(FILTER) '-' roiPrefix '-Derivative' num2str(DERIVATIVE) ...
    '-WhiteMatter' num2str(WHITE_MATTER) '-KWM' num2str(Kwm) '-CSF' num2str(CSF) '-KCSF' num2str(Kcsf) '-STANDARDIZE' num2str(STANDARDIZE) '-TSHIFT' num2str(TSHIFT) '-Echo' ECHOSTR '-' date]


%%
% first pass through data to get the mean centroid
centroid=zeros(3,nSubjects);
for s=1:nSubjects
    s
    subjIndx=subjStrs{s};
    basePath=['../data/' subjIndx '/NII/'];
    brainMaskFilename=fullfile(basePath,'resampled_brain_mask+tlrc.BRIK');
    roiMaskFilename=fullfile(basePath,['resampled_' roiPrefix '+tlrc.BRIK']);
    muRoiMaskFilename=fullfile(basePath,['resampled_' muRoiPrefix '+tlrc.BRIK']);
    [~, roiMask, ~, ~] = BrikLoad (roiMaskFilename);
    if isempty(roiMask)
        [~, roiMask, ~, ~] = BrikLoad (muRoiMaskFilename);
    end
    roiMask=logical(roiMask);
    [~, brainMask, ~, ~] = BrikLoad (brainMaskFilename);
    brainMask=logical(brainMask);
    finalMask=roiMask>0 & brainMask>0;
    [X,Y,Z]=ndgrid(1:nX,1:nY,1:nZ);
    xmask=X(finalMask); ymask=Y(finalMask); zmask=Z(finalMask);
    centroid(:,s)=[mean(xmask);mean(ymask);mean(zmask)];
end
muCentroid=mean(centroid,2); % try to align everyone's bold to for a coherent spatial average
cshifts=round(centroid-repmat(muCentroid,1,nSubjects));

%%
% main pass
alloBolds=cell(nSubjects,1);
muwBold=zeros(nX,nY,nZ,nTR);
alltShifts=zeros(nSubjects,1);
for s=1:numel(subjStrs)
    s
    subjIndx=subjStrs{s};
    basePath=['../data/' subjIndx '/NII/'];
    inputBoldFilename=fullfile(basePath,['dsbold_e' ECHOSTR '_tlrc_al+tlrc.BRIK']);
    brainMaskFilename=fullfile(basePath,'resampled_brain_mask+tlrc.BRIK');
    muRoiMaskFilename=fullfile(basePath,['resampled_' muRoiPrefix '+tlrc.BRIK']);
    roiMaskFilename=fullfile(basePath,['resampled_' roiPrefix '+tlrc.BRIK']);
    motionFilename=fullfile(basePath,['dsbold_e' ECHOSTR '_vr_motion.1D']);
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
    %[~,~,~,nTR]=size(inputBold);
    iwBoldTs=vol2ts(inputBold);
    
    [~, roiMask, ~, ~] = BrikLoad (roiMaskFilename); 
    if isempty(roiMask)
        [~, roiMask, ~, ~] = BrikLoad (muRoiMaskFilename);
    end
    roiMask=logical(roiMask);
    
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
        [~,seg,~,~] = BrikLoad (resampledSegFilename);
        wmMask=seg==WM_INDEX;
        wmBoldTs=vol2ts(outputBold,wmMask);
        [U,S,V]=svd(wmBoldTs,0);
        brainBoldTs=vol2ts(outputBold,brainMask); % only regress in the brain mask
        obrainBoldTs = regressOut(brainBoldTs.',U(:,1:Kwm).',1).';
    end
    
    if CSF
        outputBold=ts2vol(obrainBoldTs,brainMask); % put back to full 4D
        csfMask=seg==CSF_INDEX;
        csfBoldTs=vol2ts(outputBold,csfMask);
        [U,S,V]=svd(csfBoldTs,0);
        obrainBoldTs = regressOut(obrainBoldTs.',U(:,1:Kcsf).',1).';
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
    [err,ErrMessage,Info]=WriteBrikWrap(basePath,final4D,iInfo,['oBoldEcho' ECHOSTR],'tlrc');
%     
    % here we center the bold so that everyone's ROI lines up
    % this does not take into account different laser orientations
    if CENTERBOLDS
        final4Ds=circshift(final4D,[cshifts(1,s) cshifts(2,s) cshifts(3,s) 0]);
        muwBold=muwBold+final4Ds/nSubjects;
    else
        muwBold=muwBold+final4D/nSubjects;
    end
    
    % the next two steps are to get the time series of the ROI (we have
    % whole brain already)
    oBoldTs=vol2ts(final4D,finalMask);
    alloBolds{s}=oBoldTs;
end

%%
allBolds2D=cat(2,alloBolds{:});
figure;
subplot(211);
time=(0:nTR-1)*TR;
plot(time,mean(allBolds2D,2),'k');

print('-dpng',figFilenameStr);

save(outFilenameStr,'alloBolds','fl','nPad','bf','af','Kwm','Kcsf','FILTER','DERIVATIVE','WHITE_MATTER','CSF','TR','fs');

%%
% write out the average bold here
%[~, inputBold, info, ~] = BrikLoad (inputBoldFilename); % just to get "Info"
[err,ErrMessage,info]=WriteBrikWrap('../data/SAVG/NII/',muwBold,iInfo,['muwBoldEcho' ECHOSTR],'tlrc');

%%
% smooth the average bold here
for f=1:numel(FWHMSTR)
    tFWHMSTR=FWHMSTR{f};
    origPath=pwd;
    cd('../data/SAVG/NII/');
    muwPrefix=['s' tFWHMSTR 'muwBoldEcho' ECHOSTR '-WhiteMatter' num2str(WHITE_MATTER) '-KWHM' num2str(Kwm) '-CSF' num2str(CSF) '-KCSF' num2str(Kcsf) '-Center' num2str(CENTERBOLDS)];
    str=['!3dBlurInMask -prefix ' muwPrefix ' -input muwBoldEcho' ECHOSTR '+tlrc -FWHM ' tFWHMSTR ' -mask muwBoldEcho' ECHOSTR '+tlrc'];
    eval(str);
    cd(origPath);
end
% 
% %%
% % load in the smoothed average bold
% [~, smuBold, info, ~] = BrikLoad (['../data/SAVG/NII/s' FWHMSTR 'muwBoldEcho' ECHOSTR '+tlrc']);
% tmp=squeeze(smuBold(round(muCentroid(1)),round(muCentroid(2)),round(muCentroid(3)),:));
% figure; plot(time,tmp,'k');