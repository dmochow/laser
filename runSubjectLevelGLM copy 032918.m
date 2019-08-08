clear all; close all; clc
addpath(genpath('~/PROJECTS/COMMON'));

runWholeBrain=0;  % if you want to look at not just the ROI but the whole brain
useAverageROI=0; % if we want to apply the same average ROI to all subjects

% point to the data
%subjStr={'S03'};
subjStr={'S02','S03','S04','S05','S06','S07','S09','S10','S11','S12','S13'};
nSubjects=numel(subjStr);
figFilename=['../figures/' num2str(nSubjects) 'subs_grand_average ' date];

maskFilename='resampled_roi+tlrc';
controlMaskFilename='resampled_control_roi+tlrc';
altMaskFilename='muroi+tlrc'; % if subject doesn't have markers
altControlMaskFilename='muControlRoi+tlrc'; % if subject doesn't have markers
anatFilename='resampled_brain_mask+tlrc';
epiFilename_1='smnudsbold_e1_tlrc_al+tlrc';
epiFilename_2='smnudsbold_e2_tlrc_al+tlrc';
epiFilename_3='smnudsbold_e3_tlrc_al+tlrc';
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
isSigUncAll=cell(nSubjects,nEchos);
allBetas=cell(nSubjects,nEchos);
allBoldsRoi=cell(nSubjects,1);


%%
for s=1:nSubjects
    s
    path=['../data/' subjStr{s} '/NII/'];
    biopacPath=['../data/' subjStr{s} '/BIOPAC/'];
    pval=[];
    isSig_fdr=[];
    control_pval=[];

    % read in the data
    [err, mask, InfoMask, ErrMessage] = BrikLoad (fullfile(path,maskFilename));
    if err || useAverageROI
        [err, mask, InfoMask, ErrMessage] = BrikLoad (fullfile(path,altMaskFilename));
    end
    [err, controlMask, InfoMask, ErrMessage] = BrikLoad (fullfile(path,controlMaskFilename));
    if err || useAverageROI
        [err, controlMask, InfoMask, ErrMessage] = BrikLoad (fullfile(path,altControlMaskFilename));
    end
    [err, anat, Info, ErrMessage] = BrikLoad (fullfile(path,anatFilename));
    [err, bold_1, Info, ErrMessage] = BrikLoad (fullfile(path,epiFilename_1));
    [err, bold_2, Info, ErrMessage] = BrikLoad (fullfile(path,epiFilename_2));
    [err, bold_3, Info, ErrMessage] = BrikLoad (fullfile(path,epiFilename_3));
    
    % compute the mask
    if ~runWholeBrain
        finalMask=mask>0 & anat>0; % in the laser beam AND in the brain mask
        finalControlMask=controlMask>0 & anat>0; 
    else
        finalMask=anat>0;
    end
    maskInds=find(finalMask);
    controlMaskInds=find(finalControlMask);
    
    % grab the bold in the ROI
    tmp=permute(bold_1,[4 1 2 3]);
    bold(1,:,:)=tmp(:,:);
    
    tmp=permute(bold_2,[4 1 2 3]);
    bold(2,:,:)=tmp(:,:);
    
    tmp=permute(bold_3,[4 1 2 3]);
    bold(3,:,:)=tmp(:,:);
    
    bold_roi=bold(:,:,finalMask);
    allBoldsRoi{s,1}=bold_roi;
    
    bold_control_roi=bold(:,:,finalControlMask);
    allBoldsControlRoi{s,1}=bold_control_roi;
    
    %% laser onset times from biopac
    load(fullfile(biopacPath,biopacFilename),'data');
    
    [~,maxindsample]=max(diff(data(:,1)));
    onsetTimeSec=maxindsample/1000;
    onsetTimeTR=round(onsetTimeSec/TR);
    
    [~,minindsample]=min(diff(data(:,1)));
    offsetTimeSec=minindsample/1000;
    offsetTimeTR=round(offsetTimeSec/TR);
    
    if offsetTimeTR==onsetTimeTR % error
        offsetTimeTR=onsetTimeTR+214;
    end
    
    allOnsets(s,1)=onsetTimeTR;
    allOnsets(s,2)=offsetTimeTR;

    %% start GLM here (real ROI)
    % for each voxel in ROI, convert to 1D time series
    nROIvoxels=size(bold_roi,3);
    nTRs=size(bold_roi,2);
    for e=1:nEchos
        betas=[];
        for v=1:nROIvoxels
            [e,v];
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
            theseBetas=B(1:end-1); % don't include the intercept term
            betas=cat(2,betas,theseBetas);
            
            % reduced model
            Xr=ones(nTRs,1);
            Br=pinv(Xr)*Y;
            sse0=sum((Y-Xr*Br).^2);
            dof0=nTRs-size(Xr,2);
            Fstat=((sse0-sse)/(dof0-dof)) / (sse/dof);
            pval(e,v) = 1-fcdf(Fstat,size(X,2)-1,nTRs-size(X,2));
            
            %%
        end
        
        allBetas{s,e}=betas;
    end
    
    [p_fdr, p_masked] = fdr( pval, alpha);
    isSig_fdr=p_masked;
    isSig_unc=pval<alpha;
    
    isSig_e1=zeros(size(finalMask));
    isSig_e1(maskInds)=isSig_fdr(1,:);
    
    isSig_e2=zeros(size(finalMask));
    isSig_e2(maskInds)=isSig_fdr(2,:);
    
    isSig_e3=zeros(size(finalMask));
    isSig_e3(maskInds)=isSig_fdr(3,:);
    
    isSigUnc_e1=zeros(size(finalMask));
    isSigUnc_e1(maskInds)=isSig_unc(1,:);
    
    isSigUnc_e2=zeros(size(finalMask));
    isSigUnc_e2(maskInds)=isSig_unc(2,:);
    
    isSigUnc_e3=zeros(size(finalMask));
    isSigUnc_e3(maskInds)=isSig_unc(3,:);
    
    isSigAll{s,1}=isSig_e1;
    isSigAll{s,2}=isSig_e2;
    isSigAll{s,3}=isSig_e3;
    nSigVoxelsFdr=cellfun(@(x)sum(x(:)),isSigAll)
    
    isSigUncAll{s,1}=isSigUnc_e1;
    isSigUncAll{s,2}=isSigUnc_e2;
    isSigUncAll{s,3}=isSigUnc_e3;
    nSigVoxelsUnc=cellfun(@(x)sum(x(:)),isSigUncAll)
    
    if exist(fullfile(path,'isSig_e1+tlrc.HEAD'),'file')==2
        delete(fullfile(path,'isSig_e1+tlrc.HEAD'));
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
    
    
    if exist(fullfile(path,'isSig_e2+tlrc.HEAD'),'file')==2
        delete(fullfile(path,'isSig_e2+tlrc.HEAD'));
        delete(fullfile(path,'isSig_e2+tlrc.BRIK'));
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
    
    if exist(fullfile(path,'isSig_e3+tlrc.HEAD'),'file')==2
        delete(fullfile(path,'isSig_e3+tlrc.HEAD'));
        delete(fullfile(path,'isSig_e3+tlrc.BRIK'));
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
    
    %% control GLM starts here 
    nControlROIvoxels=size(bold_control_roi,3);
    for e=1:nEchos
        controlBetas=[];
        for v=1:nControlROIvoxels
            [e,v];
            Y=squeeze(bold_control_roi(e,:,v));
            Y=Y(:);
            
            % construct design matrix
            %%
            X=eye(nTRs);
            X=X(:,onsetTimeTR:offsetTimeTR);
            X=cat(2,X,ones(nTRs,1));
            B=pinv(X)*Y;
            sse=sum((Y-X*B).^2);
            dof=nTRs-size(X,2);

            theseControlBetas=B(1:end-1);
            controlBetas=cat(2,controlBetas,theseControlBetas);
            
            % reduced model
            Xr=ones(nTRs,1);
            Br=pinv(Xr)*Y;
            sse0=sum((Y-Xr*Br).^2);
            dof0=nTRs-size(Xr,2);
            Fstat=((sse0-sse)/(dof0-dof)) / (sse/dof);
            control_pval(e,v) = 1-fcdf(Fstat,size(X,2)-1,nTRs-size(X,2));
            %isSig_bonf(e,v)=pval(e,v)<alpha/prod(size(pval));
            %%
        end
        
        allControlBetas{s,e}=controlBetas;
    end
    
    [p_fdr_control, p_masked_control] = fdr( control_pval, alpha);
    isSig_fdr_control=p_masked_control;
    isSig_unc_control=control_pval<alpha;
    
    isSig_e1_control=zeros(size(finalControlMask));
    isSig_e1_control(controlMaskInds)=isSig_fdr_control(1,:);
    
    isSig_e2_control=zeros(size(finalControlMask));
    isSig_e2_control(controlMaskInds)=isSig_fdr_control(2,:);
    
    isSig_e3_control=zeros(size(finalControlMask));
    isSig_e3_control(controlMaskInds)=isSig_fdr_control(3,:);
    
    isSigUnc_e1_control=zeros(size(finalControlMask));
    isSigUnc_e1_control(controlMaskInds)=isSig_unc_control(1,:);
    
    isSigUnc_e2_control=zeros(size(finalControlMask));
    isSigUnc_e2_control(controlMaskInds)=isSig_unc_control(2,:);
    
    isSigUnc_e3_control=zeros(size(finalControlMask));
    isSigUnc_e3_control(controlMaskInds)=isSig_unc_control(3,:);
    
    isSigAllControl{s,1}=isSig_e1_control;
    isSigAllControl{s,2}=isSig_e2_control;
    isSigAllControl{s,3}=isSig_e3_control;
    nSigVoxelsFdrControl=cellfun(@(x)sum(x(:)),isSigAllControl)
    
    isSigUncAllControl{s,1}=isSigUnc_e1_control;
    isSigUncAllControl{s,2}=isSigUnc_e2_control;
    isSigUncAllControl{s,3}=isSigUnc_e3_control;
    nSigVoxelsUncControl=cellfun(@(x)sum(x(:)),isSigUncAllControl)
    
    if exist(fullfile(path,'isSig_e1_control+tlrc.HEAD'),'file')==2
        delete(fullfile(path,'isSig_e1_control+tlrc.HEAD'));
        delete(fullfile(path,'isSig_e1_control+tlrc.BRIK'));
    end
    
    InfoWrite=InfoMask;
    opt.Scale=0;
    opt.Prefix='isSig_e1_control';
    opt.View='tlrc';
    opt.Verbose=1;
    opt.AppendHistory=1;
    opt.NoCheck=0;
    opt.Overwrite=1;
    origPath=pwd;
    cd(path);
    [err, ErrMessage, Info] = WriteBrik (isSig_e1_control, InfoWrite,opt);
    cd(origPath);
    
    
    if exist(fullfile(path,'isSig_e2_control+tlrc.HEAD'),'file')==2
        delete(fullfile(path,'isSig_e2_control+tlrc.HEAD'));
        delete(fullfile(path,'isSig_e2_control+tlrc.BRIK'));
    end
    
    
    InfoWrite=InfoMask;
    opt.Scale=0;
    opt.Prefix='isSig_e2_control';
    opt.View='tlrc';
    opt.Verbose=1;
    opt.AppendHistory=1;
    opt.NoCheck=0;
    opt.Overwrite=1;
    origPath=pwd;
    cd(path);
    [err, ErrMessage, Info] = WriteBrik (isSig_e2_control, InfoWrite,opt);
    cd(origPath);
    
    if exist(fullfile(path,'isSig_e3_control+tlrc.HEAD'),'file')==2
        delete(fullfile(path,'isSig_e3_control+tlrc.HEAD'));
        delete(fullfile(path,'isSig_e3_control+tlrc.BRIK'));
    end
    
    
    InfoWrite=InfoMask;
    opt.Scale=0;
    opt.Prefix='isSig_e3_control';
    opt.View='tlrc';
    opt.Verbose=1;
    opt.AppendHistory=1;
    opt.NoCheck=0;
    opt.Overwrite=1;
    origPath=pwd;
    cd(path);
    [err, ErrMessage, Info] = WriteBrik (isSig_e3_control, InfoWrite,opt);
    cd(origPath);
    
end



%%
% different subjects have different number of TRs during laser on, so
% figure out the minimum time that laser was on for ALL subjects
for s=1:nSubjects
    for e=1:nEchos
        nTRs_laser(s,e)=size(allBetas{s,e},1);
    end
end
minTRs=min(nTRs_laser(:));
catBetas1=[];
catBetas2=[];
catBetas3=[];
for s=1:nSubjects
%         catBetas1=cat(2,catBetas1,allBetas{s,1}(1:minTRs,:));
%         catBetas2=cat(2,catBetas2,allBetas{s,2}(1:minTRs,:));
%         catBetas3=cat(2,catBetas3,allBetas{s,3}(1:minTRs,:));
        catBetas1=cat(2,catBetas1,allControlBetas{s,1}(1:minTRs,:));
        catBetas2=cat(2,catBetas2,allControlBetas{s,2}(1:minTRs,:));
        catBetas3=cat(2,catBetas3,allControlBetas{s,3}(1:minTRs,:));
end
muBetas(:,1)=mean(catBetas1,2);
muBetas(:,2)=mean(catBetas2,2);
muBetas(:,3)=mean(catBetas3,2);

figure;
for e=1:3
    subplot(nEchos,1,e);
    plot((0:minTRs-1)*TR,muBetas(:,e));
    %xlim([0 1800]);
end
%%
% because I'm eager
time=(0:644)*TR;
tLaser=[600 1200];
grandBold=mean(cat(3,allBoldsRoi{:}),3);
figure; 
subplot(411); hold on
plot(time,grandBold(1,:)); title('Echo 1');
yl=ylim;
harea=area(tLaser,yl(2)*[1 1],'BaseValue',yl(1));
set(harea,'FaceColor',[0.7 0.7 0.7]);
set(harea,'FaceAlpha',0.5);
set(harea,'EdgeColor','none');
subplot(412); hold on
plot(time,grandBold(2,:)); title('Echo 2');
yl=ylim;
harea=area(tLaser,yl(2)*[1 1],'BaseValue',yl(1));
set(harea,'FaceColor',[0.7 0.7 0.7]);
set(harea,'FaceAlpha',0.5);
set(harea,'EdgeColor','none');
subplot(413); hold on
plot(time,grandBold(3,:)); title('Echo 3');
yl=ylim;
harea=area(tLaser,yl(2)*[1 1],'BaseValue',yl(1));
set(harea,'FaceColor',[0.7 0.7 0.7]);
set(harea,'FaceAlpha',0.5);
set(harea,'EdgeColor','none');
subplot(414); hold on
plot(time,mean(grandBold)); title('Echo mean');
yl=ylim;
harea=area(tLaser,yl(2)*[1 1],'BaseValue',yl(1));
set(harea,'FaceColor',[0.7 0.7 0.7]);
set(harea,'FaceAlpha',0.5);
set(harea,'EdgeColor','none');
print('-depsc','-r600',figFilename);

%%
% because I'm eager
time=(0:644)*TR;
tLaser=[600 1200];
grandBold=mean(cat(3,allBoldsControlRoi{:}),3);
figure; 
subplot(411); hold on
plot(time,grandBold(1,:)); title('Echo 1');
yl=ylim;
harea=area(tLaser,yl(2)*[1 1],'BaseValue',yl(1));
set(harea,'FaceColor',[0.7 0.7 0.7]);
set(harea,'FaceAlpha',0.5);
set(harea,'EdgeColor','none');
subplot(412); hold on
plot(time,grandBold(2,:)); title('Echo 2');
yl=ylim;
harea=area(tLaser,yl(2)*[1 1],'BaseValue',yl(1));
set(harea,'FaceColor',[0.7 0.7 0.7]);
set(harea,'FaceAlpha',0.5);
set(harea,'EdgeColor','none');
subplot(413); hold on
plot(time,grandBold(3,:)); title('Echo 3');
yl=ylim;
harea=area(tLaser,yl(2)*[1 1],'BaseValue',yl(1));
set(harea,'FaceColor',[0.7 0.7 0.7]);
set(harea,'FaceAlpha',0.5);
set(harea,'EdgeColor','none');
subplot(414); hold on
plot(time,mean(grandBold)); title('Echo mean');
yl=ylim;
harea=area(tLaser,yl(2)*[1 1],'BaseValue',yl(1));
set(harea,'FaceColor',[0.7 0.7 0.7]);
set(harea,'FaceAlpha',0.5);
set(harea,'EdgeColor','none');
%print('-depsc','-r600',figFilename);
%%
save(['../data/precomputed/subjectLevelGlmResults ' date],'nSigVoxelsFdr','nSigVoxelsUnc','allBoldsRoi','allBetas','allOnsets','nSigVoxelsFdr','nSigVoxelsUnc','useAverageROI','runWholeBrain','subjStr');

