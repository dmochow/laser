% 10.17.18 estimating T2s from raw bold (dsbold_) to get exponential shape
% this has to be done subject wise
% 10.16.18
% optimal multi echo combination following Kundu 2012 following Posse 1999
%
clear all; close all; clc
subjStrs={'S02','S03','S04','S05','S06','S07','S09','S10','S11','S12','S13'...
    ,'S14','S15','S16','S17','S18','S19','S20','S21','S22'};
nSubjects=numel(subjStrs);
nEchos=3;
nTRs=645;
TEs=[12.8,34.3,55.6]'/1000;

for s=1:nSubjects
    fprintf('Processing subject %d of %d \n',s,nSubjects);
    
    pathToData=['../data/' subjStrs{s} '/NII'];
    
    %
    brainMaskFilename='resampled_brain_mask+tlrc';
    [~, brainMask, info, ~] = BrikLoad (fullfile(pathToData,brainMaskFilename));
    brainMask=logical(brainMask);
    
    % preprocessed but not denoised bold
    for e=1:nEchos
        ECHOSTR=num2str(e);
        inputBoldFilename=fullfile(pathToData,['dsbold_e' ECHOSTR '_tlrc_al+tlrc.BRIK']);
        [~, tBold, infoBold, ~] = BrikLoad (inputBoldFilename);
        boldTs(:,:,e)=vol2ts(tBold,brainMask);
    end
    muBoldTs=squeeze(mean(boldTs,1));
    logMuBoldTs=log(muBoldTs); % log-linear regression
    
    %
    nVoxels=size(muBoldTs,1);
    R2=zeros(nVoxels,1);
    logSo=zeros(nVoxels,1);
    for v=1:nVoxels
        v
        [b,bint,r,rint,stats] = regress(logMuBoldTs(v,:).',[TEs ones(3,1)]);
        logSo(v)=b(2);
        R2(v,1)=-b(1);
    end
    
    % fix infs in R2 estimate
    muR2=nanmean(R2(~isinf(R2)));
    R2(isinf(R2))=muR2;
    
    %R2c=R2;
    %R2c(isinf(R2))=NaN;
    
    %%
    RR2=repmat(R2,1,nEchos);
    TTE=repmat(TEs.',nVoxels,1);
    W=TTE.*exp(-TTE.*RR2);
    W=W./repmat(sum(W,2),1,nEchos);
    
    %%
    W3=zeros(nTRs,size(W,1),size(W,2));
    for t=1:nTRs
        W3(t,:,:)=W;
    end
    
    %%
    oef=sum(boldTs.*W3,3);
    
    oef4D=ts2vol(oef,brainMask);
    
    % write the OEF weighted BOLD
    [err,ErrMessage,info]=WriteBrikWrap(pathToData,oef4D,infoBold,'oc_dsbold','tlrc');
    
    clear boldTs; % changes size with subject
    
end
