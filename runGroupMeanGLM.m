clear all; close all; clc
addpath(genpath('~/PROJECTS/COMMON'));

% point to the data
subjStr={'SAVG'};
nSubjects=numel(subjStr);
figFilename=['../figures/subject_average_grand_average'];

maskFilename='resampled_roi+tlrc';
altMaskFilename='muroi+tlrc';
anatFilename='resampled_TT_N27+tlrc';
epiFilename_1='groupMean_e1+tlrc.BRIK';
epiFilename_2='groupMean_e2+tlrc.BRIK';
epiFilename_3='groupMean_e3+tlrc.BRIK';
biopacFilename='biopac.mat';
TEs=[12.8 34.13 55.46];
nEchos=numel(TEs);
TR=2.8;
alpha=0.05;

allBolds=cell(nSubjects,1);
allBoldsControl=cell(nSubjects,1);
allMaskInds=cell(nSubjects,1);
allOnsets=zeros(nSubjects,2);

isSigAll=cell(nSubjects,nEchos);

for s=1:nSubjects
    s
    path=['../data/' subjStr{s} '/'];
    biopacPath=['../data/' subjStr{s} '/BIOPAC/'];
    pval=[];
    isSig_fdr=[];
    isSig_fdr=[];
    
    % read in the data
    [err, mask, InfoMask, ErrMessage] = BrikLoad (fullfile(path,altMaskFilename));
    [err, anat, Info, ErrMessage] = BrikLoad (fullfile(path,anatFilename));
    [err, bold_1, Info, ErrMessage] = BrikLoad (fullfile(path,epiFilename_1));
    [err, bold_2, Info, ErrMessage] = BrikLoad (fullfile(path,epiFilename_2));
    [err, bold_3, Info, ErrMessage] = BrikLoad (fullfile(path,epiFilename_3));
    
    % compute the mask
    finalMask=mask>0 & anat>0; % in the laser beam AND in the brain
    maskInds=find(finalMask);
    
    % grab the bold in the ROI
    tmp=permute(bold_1,[4 1 2 3]);
    bold(1,:,:)=tmp(:,:);
    
    tmp=permute(bold_2,[4 1 2 3]);
    bold(2,:,:)=tmp(:,:);
    
    tmp=permute(bold_3,[4 1 2 3]);
    bold(3,:,:)=tmp(:,:);
    
    bold_roi=bold(:,:,finalMask);
    
    % onset times
%     load(fullfile(biopacPath,biopacFilename),'data');
%     
%     [~,maxindsample]=max(diff(data(:,1)));
%     onsetTimeSec=maxindsample/1000;
%     onsetTimeTR=round(onsetTimeSec/TR);
%     
%     [~,minindsample]=min(diff(data(:,1)));
%     offsetTimeSec=minindsample/1000;
%     offsetTimeTR=round(offsetTimeSec/TR);
%     
%     if offsetTimeTR==onsetTimeTR % error
%         offsetTimeTR=onsetTimeTR+214;
%     end
%     
%     allOnsets(s,1)=onsetTimeTR;
%     allOnsets(s,2)=offsetTimeTR;

    onsetTimeTR=214;
    offsetTimeTR=429;
    
    % for each voxel in ROI, convert to 1D time series
    nROIvoxels=size(bold_roi,3);
    nTRs=size(bold_roi,2);
    for e=1:nEchos
        for v=1:nROIvoxels
            [e,v]
            Y=squeeze(bold_roi(e,:,v));
            Y=Y(:);
            
            % construct design matrix
            %%
            X=eye(nTRs);
            X=X(:,onsetTimeTR:offsetTimeTR);
            X=cat(2,X,ones(nTRs,1));
            % Y=X*B --> B = pinv(X)*Y
            B=pinv(X)*Y;
            sse=sum((Y-X*B).^2);
            dof=nTRs-size(X,2);
            
            % reduced model
            Xr=ones(nTRs,1);
            Br=pinv(Xr)*Y;
            sse0=sum(Y.^2);
            dof0=nTRs-size(Xr,2);
            Fstat=((sse0-sse)/(dof0-dof)) / (sse/dof);
            pval(e,v) = 1-fcdf(Fstat,size(X,2)-1,nTRs-size(X,2));
            %isSig_bonf(e,v)=pval(e,v)<alpha/prod(size(pval));            
            
            
            %%
        end
    end
    
    [p_fdr, p_masked] = fdr( pval, alpha);
    isSig_fdr=p_masked;
    
    isSig_e1=zeros(size(finalMask));
    isSig_e1(maskInds)=isSig_fdr(1,:);
    
    isSig_e2=zeros(size(finalMask));
    isSig_e2(maskInds)=isSig_fdr(2,:);
    
    isSig_e3=zeros(size(finalMask));
    isSig_e3(maskInds)=isSig_fdr(3,:);
    
    isSigAll{s,1}=isSig_e1;
    isSigAll{s,2}=isSig_e2;
    isSigAll{s,3}=isSig_e3;
    tmp=cellfun(@(x)sum(x(:)),isSigAll)
    
    if exist(fullfile(path,'isSig_e1+tlrc.BRIK'),'file')
        %delete isSig_e1+tlrc.BRIK;
        delete(fullfile(path,'isSig_e1+tlrc.BRIK'));
        %delete isSig_e1+tlrc.HEAD;
        delete(fullfile(path,'isSig_e1+tlrc.BRIK'));
    end
    
    InfoWrite=InfoMask;
    opt.Scale=0;
    opt.Prefix='isSig_e1';
    opt.View='tlrc';
    opt.Verbose=1;
    opt.AppendHistory=1;
    opt.NoCheck=0;
    opt.Overwrite=1;
    origPath=pwd;
    cd(path);
    [err, ErrMessage, Info] = WriteBrik (isSig_e1, InfoWrite,opt);
    cd(origPath);
    
    
    if exist(fullfile(path,'isSig_e2+tlrc.BRIK'),'file')
        delete isSig_e2+tlrc.BRIK;
        delete isSig_e2+tlrc.HEAD;
    end
    InfoWrite=InfoMask;
    opt.Scale=0;
    opt.Prefix='isSig_e2';
    opt.View='tlrc';
    opt.Verbose=1;
    opt.AppendHistory=1;
    opt.NoCheck=0;
    opt.Overwrite=1;
    origPath=pwd;
    cd(path);
    [err, ErrMessage, Info] = WriteBrik (isSig_e2, InfoWrite,opt);
    cd(origPath);
    

    if exist(fullfile(path,'isSig_e3+tlrc.BRIK'),'file')
        delete isSig_e3+tlrc.BRIK;
        delete isSig_e3+tlrc.HEAD;
    end
    InfoWrite=InfoMask;
    opt.Scale=0;
    opt.Prefix='isSig_e3';
    opt.View='tlrc';
    opt.Verbose=1;
    opt.AppendHistory=1;
    opt.NoCheck=0;
    opt.Overwrite=1;
    origPath=pwd;
    cd(path);
    [err, ErrMessage, Info] = WriteBrik (isSig_e3, InfoWrite,opt);
    cd(origPath);
    
end

%%


