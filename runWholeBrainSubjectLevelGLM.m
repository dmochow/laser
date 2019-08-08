clear all; close all; clc
addpath(genpath('~/PROJECTS/COMMON'));

useAverageROI=0; % if we want to apply the same average ROI to all subjects
% point to the data
subjStr={'S03'};
%subjStr={'S02','S03','S04','S05','S06','S07','S09','S10','S11','S12','S13'};
nSubjects=numel(subjStr);

brainMaskFilename='resampled_brain_mask+tlrc';
epiFilename_1='smnudsbold_e1_tlrc_al+tlrc';
epiFilename_2='smnudsbold_e2_tlrc_al+tlrc';
epiFilename_3='smnudsbold_e3_tlrc_al+tlrc';
biopacFilename='biopac.mat';
TEs=[12.8 34.13 55.46];
nEchos=numel(TEs);
TR=2.8;
alpha=0.05;

allBolds=cell(nSubjects,1);
allOnsets=zeros(nSubjects,2);
isSigAll=cell(nSubjects,nEchos);
isSigUncAll=cell(nSubjects,nEchos);
allBetas=cell(nSubjects,1);

%%
for s=1:nSubjects
    s
    path=['../data/' subjStr{s} '/NII/'];
    biopacPath=['../data/' subjStr{s} '/BIOPAC/'];

    % read in the data
    [err, brainMask, Info, ErrMessage] = BrikLoad (fullfile(path,brainMaskFilename));
    [err, bold_1, Info, ErrMessage] = BrikLoad (fullfile(path,epiFilename_1));
    [err, bold_2, Info, ErrMessage] = BrikLoad (fullfile(path,epiFilename_2));
    [err, bold_3, Info, ErrMessage] = BrikLoad (fullfile(path,epiFilename_3));
    
    % compute the mask
    finalMask=brainMask>0;

    % grab the bold in the ROI
    tmp=permute(bold_1,[4 1 2 3]);
    bold(1,:,:)=tmp(:,:);
    
    tmp=permute(bold_2,[4 1 2 3]);
    bold(2,:,:)=tmp(:,:);
    
    tmp=permute(bold_3,[4 1 2 3]);
    bold(3,:,:)=tmp(:,:);
    
    boldMasked=bold(:,:,finalMask);
    [~,nTRs,nVoxels]=size(boldMasked);
    
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
    nTRsLaser=offsetTimeTR-onsetTimeTR+1;

    % start GLM here
    pval=zeros(nEchos,nVoxels);
    betas=zeros(nEchos,nVoxels,nTRsLaser);
    parfor e=1:nEchos
        for v=1:nVoxels
            [e,v/nVoxels*100]
            Y=squeeze(boldMasked(e,:,v));
            Y=Y(:);
            Xf=eye(nTRs);
            Xf=Xf(:,onsetTimeTR:offsetTimeTR);
            statsf=myGLM(Y,Xf);
            ssef=statsf.sse;
            doff=statsf.dof;
            
            betas(e,v,:)=statsf.B(1:end-1);
            
            % reduced model
            Xr=ones(nTRs,1);
            statsr=myGLM(Y,Xr);
            sser=statsr.sse;
            dofr=statsr.dof;
            Fstat=((sser-ssef)/(dofr-doff)) / (ssef/doff);
            pval(e,v)=fcdf(Fstat,dofr-doff,doff,'upper');
            
        end
        
        
    end
    
    allBetas{s,1}=betas;
    
    [p_fdr, p_masked] = fdr( pval, alpha);
    isSig_fdr=p_masked;
    
    isSig_e1=zeros(size(finalMask));
    isSig_e1(finalMask)=isSig_fdr(1,:);
    
    isSig_e2=zeros(size(finalMask));
    isSig_e2(finalMask)=isSig_fdr(2,:);
    
    isSig_e3=zeros(size(finalMask));
    isSig_e3(finalMask)=isSig_fdr(3,:);
    
    isSigAll{s,1}=isSig_e1;
    isSigAll{s,2}=isSig_e2;
    isSigAll{s,3}=isSig_e3;
    nSigVoxelsFdr=cellfun(@(x)sum(x(:)),isSigAll)
    
    if exist(fullfile(path,'isSig_e1+tlrc.HEAD'),'file')==2
        delete(fullfile(path,'isSig_e1+tlrc.HEAD'));
        delete(fullfile(path,'isSig_e1+tlrc.BRIK'));
    end
    
    if exist(fullfile(path,'isSig_e2+tlrc.HEAD'),'file')==2
        delete(fullfile(path,'isSig_e2+tlrc.HEAD'));
        delete(fullfile(path,'isSig_e2+tlrc.BRIK'));
    end
      
    if exist(fullfile(path,'isSig_e3+tlrc.HEAD'),'file')==2
        delete(fullfile(path,'isSig_e3+tlrc.HEAD'));
        delete(fullfile(path,'isSig_e3+tlrc.BRIK'));
    end
    
    WriteBrikWrap(path,isSig_e1,Info,'isSig_e1','tlrc');
    WriteBrikWrap(path,isSig_e2,Info,'isSig_e2','tlrc');
    WriteBrikWrap(path,isSig_e3,Info,'isSig_e3','tlrc');
    
end



% %%
% % different subjects have different number of TRs during laser on, so
% % figure out the minimum time that laser was on for ALL subjects
% for s=1:nSubjects
%     for e=1:nEchos
%         nTRs_laser(s,e)=size(allBetas{s,e},1);
%     end
% end
% minTRs=min(nTRs_laser(:));
% catBetas1=[];
% catBetas2=[];
% catBetas3=[];
% for s=1:nSubjects
% %         catBetas1=cat(2,catBetas1,allBetas{s,1}(1:minTRs,:));
% %         catBetas2=cat(2,catBetas2,allBetas{s,2}(1:minTRs,:));
% %         catBetas3=cat(2,catBetas3,allBetas{s,3}(1:minTRs,:));
%         catBetas1=cat(2,catBetas1,allControlBetas{s,1}(1:minTRs,:));
%         catBetas2=cat(2,catBetas2,allControlBetas{s,2}(1:minTRs,:));
%         catBetas3=cat(2,catBetas3,allControlBetas{s,3}(1:minTRs,:));
% end
% muBetas(:,1)=mean(catBetas1,2);
% muBetas(:,2)=mean(catBetas2,2);
% muBetas(:,3)=mean(catBetas3,2);
% 
% figure;
% for e=1:3
%     subplot(nEchos,1,e);
%     plot((0:minTRs-1)*TR,muBetas(:,e));
%     %xlim([0 1800]);
% end
% %%
% % because I'm eager
% time=(0:644)*TR;
% tLaser=[600 1200];
% grandBold=mean(cat(3,allBoldsRoi{:}),3);
% figure; 
% subplot(411); hold on
% plot(time,grandBold(1,:)); title('Echo 1');
% yl=ylim;
% harea=area(tLaser,yl(2)*[1 1],'BaseValue',yl(1));
% set(harea,'FaceColor',[0.7 0.7 0.7]);
% set(harea,'FaceAlpha',0.5);
% set(harea,'EdgeColor','none');
% subplot(412); hold on
% plot(time,grandBold(2,:)); title('Echo 2');
% yl=ylim;
% harea=area(tLaser,yl(2)*[1 1],'BaseValue',yl(1));
% set(harea,'FaceColor',[0.7 0.7 0.7]);
% set(harea,'FaceAlpha',0.5);
% set(harea,'EdgeColor','none');
% subplot(413); hold on
% plot(time,grandBold(3,:)); title('Echo 3');
% yl=ylim;
% harea=area(tLaser,yl(2)*[1 1],'BaseValue',yl(1));
% set(harea,'FaceColor',[0.7 0.7 0.7]);
% set(harea,'FaceAlpha',0.5);
% set(harea,'EdgeColor','none');
% subplot(414); hold on
% plot(time,mean(grandBold)); title('Echo mean');
% yl=ylim;
% harea=area(tLaser,yl(2)*[1 1],'BaseValue',yl(1));
% set(harea,'FaceColor',[0.7 0.7 0.7]);
% set(harea,'FaceAlpha',0.5);
% set(harea,'EdgeColor','none');
% print('-depsc','-r600',figFilename);
% 
% %%
% % because I'm eager
% time=(0:644)*TR;
% tLaser=[600 1200];
% grandBold=mean(cat(3,allBoldsControlRoi{:}),3);
% figure; 
% subplot(411); hold on
% plot(time,grandBold(1,:)); title('Echo 1');
% yl=ylim;
% harea=area(tLaser,yl(2)*[1 1],'BaseValue',yl(1));
% set(harea,'FaceColor',[0.7 0.7 0.7]);
% set(harea,'FaceAlpha',0.5);
% set(harea,'EdgeColor','none');
% subplot(412); hold on
% plot(time,grandBold(2,:)); title('Echo 2');
% yl=ylim;
% harea=area(tLaser,yl(2)*[1 1],'BaseValue',yl(1));
% set(harea,'FaceColor',[0.7 0.7 0.7]);
% set(harea,'FaceAlpha',0.5);
% set(harea,'EdgeColor','none');
% subplot(413); hold on
% plot(time,grandBold(3,:)); title('Echo 3');
% yl=ylim;
% harea=area(tLaser,yl(2)*[1 1],'BaseValue',yl(1));
% set(harea,'FaceColor',[0.7 0.7 0.7]);
% set(harea,'FaceAlpha',0.5);
% set(harea,'EdgeColor','none');
% subplot(414); hold on
% plot(time,mean(grandBold)); title('Echo mean');
% yl=ylim;
% harea=area(tLaser,yl(2)*[1 1],'BaseValue',yl(1));
% set(harea,'FaceColor',[0.7 0.7 0.7]);
% set(harea,'FaceAlpha',0.5);
% set(harea,'EdgeColor','none');
% %print('-depsc','-r600',figFilename);
% %%
% save(['../data/precomputed/subjectLevelGlmResults ' date],'nSigVoxelsFdr','nSigVoxelsUnc','allBoldsRoi','allBetas','allOnsets','nSigVoxelsFdr','nSigVoxelsUnc','useAverageROI','runWholeBrain','subjStr');
% 
