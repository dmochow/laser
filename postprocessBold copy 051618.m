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

subjStrs={'S02','S03','S04','S05','S06','S07','S09','S10','S11','S12','S13','S14','S15','S16','S17','S18','S19'};
nSubjects=numel(subjStrs);
TR=2.8;
fs=1/TR;
fl=0.1;
nPad=30;
[bf,af]=butter(3,fl/(fs/2));
WM=3; % segmentation label for white matter
Kwm=5; % number of WM PCs to remove

% OPTIONS
FILTER=1;   % lowpass filter at fl Hz
DERIVATIVE=1; % add the derivative of the realignment parameters to regressors
WHITE_MATTER=1; % whether to regress out the white matter PCs


alloBolds=cell(nSubjects,1);
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
    [nx,ny,nz,nTR]=size(inputBold);
     
    [~, roiMask, ~, ~] = BrikLoad (roiMaskFilename);
    roiMask=logical(roiMask);
    
    if isempty(roiMask)
        [~, roiMask, ~, ~] = BrikLoad (muRoiMaskFilename);
        roiMask=logical(roiMask);
    end
    
    [~, brainMask, ~, ~] = BrikLoad (brainMaskFilename);
    brainMask=logical(brainMask);
    
    finalMask=roiMask>0 & brainMask>0;
    
    iwBoldTs=vol2ts(inputBold);
    
    fid=fopen(motionFilename);
    mals=textscan(fid,'%f %f %f %f %f %f');
    mals=cell2mat(mals);
    if DERIVATIVE
        dmals=cat(1,zeros(1,6),diff(mals));
        mals=cat(2,mals,dmals);
    end
    mals=mals-repmat(mean(mals),size(mals,1),1);
    
    owBoldTs = regressOut(iwBoldTs.',mals.',1).';
    outputBold=ts2vol(owBoldTs,logical(ones(nx,ny,nz)));
    %myoBoldTs=vol2ts(myOutputBold,finalMask);
    
    
    if WHITE_MATTER
        % processing myOutputBold
        [~,seg,~,~] = BrikLoad (resampledSegFilename);
        wmMask=seg==WM;
        wmBoldTs=vol2ts(outputBold,wmMask);
        [U,S,V]=svd(wmBoldTs,0);
        owBoldTs = regressOut(owBoldTs.',U(:,1:Kwm).',1).';
        % regress out first K PCs
    end
    
    tmp4D=ts2vol(owBoldTs,logical(ones(nx,ny,nz)));
    oBoldTs=vol2ts(tmp4D,finalMask);
    
    if FILTER
        tmp=cat(1,repmat(oBoldTs(1,:),nPad,1),oBoldTs);
        tmpOut=filter(bf,af,tmp,[],1);
        foBoldTs=tmpOut(nPad+1:end,:);
    end
    
    alloBolds{s}=foBoldTs;
    
end

%%
allBolds2D=cat(2,alloBolds{:});
figure;
subplot(211);
time=(0:nTR-1)*TR;
plot(time,mean(allBolds2D,2),'k');

print -dpng ../figures/tmp

% %% ADD DERIVATIVES AND ABSOLUTE VALUES TO REGRESSION
% iwBoldTs=vol2ts(inputBold);
% oBoldTs=vol2ts(outputBold,roiMask);
%
% fid=fopen(motionFilename);
% mals=textscan(fid,'%f %f %f %f %f %f');
% mals=cell2mat(mals);
% dmals=cat(1,zeros(1,6),diff(mals));
% amals=abs(mals);
% mals=cat(2,cat(2,mals,dmals),amals);
% mals_z=zscore(mals);
%
% owBoldTs = regressOut(iwBoldTs.',mals_z.',1).';
% myOutputBold=ts2vol(owBoldTs,logical(ones(nx,ny,nz)));
% myoBoldTs=vol2ts(myOutputBold,roiMask);
%
% %% ADD FILTER AT 0.1 Hz
% fs=1/TR;
% [bf,af]=butter(4,0.1/(fs/2));
%
% nPad=30;
% tmp=cat(1,repmat(myoBoldTs(1,:),nPad,1),myoBoldTs);
% tmpOut=filter(bf,af,tmp,[],1);
% myfoBoldTs=tmpOut(nPad+1:end,:);
%
% figure;
% subplot(211);
% plot(mean(oBoldTs,2));
% subplot(212);
% plot(mean(myfoBoldTs,2));
% corrcoef(mean(oBoldTs,2),mean(myfoBoldTs,2))

%
% %%
% % segment skull stripped anatomical
% str=['!3dSeg -anat ' basePath 'anat_ns+orig -mask AUTO -classes ''CSF ; GM ; WM'' -bias_classes ''GM ; WM'' -bias_fwhm 25 -mixfrac UNI -main_N 5 -blur_meth BFT -prefix ' basePath];
% eval(str);
%
% %%
% % coregister segmentation with talairach
% origPath=pwd;
% cd(basePath);
% str='!@auto_tlrc -apar anat+tlrc -input Classes+orig.BRIK';
% eval(str);
% cd(origPath);
%
% %%
% % resample talairach segmentation on bold grid
% origPath=pwd;
% cd(basePath);
% str='!3dresample -master dsbold_e2_tlrc_al+tlrc -prefix resampled_Classes -input Classes+tlrc';
% eval(str);
% cd(origPath);
%
%
%
%
%
% %%
% [~, outputBold, Info, ~] = BrikLoad (outputBrik);
%
% %%

%
% % % draw

%
%
% %%
% % regress out CSF
% [~, segBold, Info, ~] = BrikLoad (segBrik);
%
% CSF=1;
% csfMask=segBold==CSF;
% csfBoldTs=vol2ts(outputBold,csfMask);
% [Uc,Sc,Vc]=svd(csfBoldTs,0);
%
% %%
% WM=3;
% wmMask=segBold==WM;
% wmBoldTs=vol2ts(outputBold,wmMask);
% [Uw,Sw,Vw]=svd(wmBoldTs,0);
%
% %%
% wbBoldTs=vol2ts(outputBold);
% ooBoldTs=regressOut(oBoldTs.',Uc(:,1:5).',1).';
% oooBoldTs=regressOut(oBoldTs.',Uw(:,1:5).',1).';
% %%
% %
%
% %% draw
% figure;
% subplot(211);
% plot(mean(oBoldTs,2));
% subplot(212);
% %plot(mean(ooBoldTs,2));
% plot(mean(oooBoldTs,2));
%
% %%
% figure; hold on
% plot(mean(csfBoldTs,2));
% plot(mean(wmBoldTs,2));
%
%
% %precomputedFilename=['../data/precomputed/allBoldsRoi 22-Apr-2018'];
% %load(precomputedFilename,'allBoldsRoiOut','allOnsets','TR','TEs','subjStr','nSubjects','nEchos');
% % sBold=mean(squeeze(allBoldsRoiOut{sIndx}(1,:,:)),2);
% % % pull in motion alignment parameters
% % malFile='../data/S02/NII/dsbold_e1_vr_motion.1D';
% % fid=fopen(malFile);
% % mals=textscan(fid,'%f %f %f %f %f %f');
% % mals=cell2mat(mals);
% % dM=cat(1,zeros(1,6),diff(mals));
% % figure
% % subplot(211);
% % plot(sBold);
% % subplot(212);
% % plot(sum(dM.^2,2));