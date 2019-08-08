% make mock bold series using autoregressive modeling
%
% examine scaling of all BOLD series (output of 3dTproject)
clear all; close all; clc
addpath(genpath('~/PROJECTS/COMMON'));

useAverageROI=0; % if we want to apply the same average ROI to all subjects
subjStr={'S02','S03','S04','S05','S06','S07','S09','S10','S11','S12','S13','S14','S15'};
nSubjects=numel(subjStr);
writeBriks=0;

maskFilename='resampled_roi+tlrc';
controlMaskFilename='resampled_control_roi+tlrc';
altMaskFilename='muroi+tlrc'; % if subject doesn't have markers
altControlMaskFilename='muControlRoi+tlrc'; % if subject doesn't have markers
brainMaskFilename='resampled_brain_mask+tlrc';
epiFilename_1='smnudsbold_e1_tlrc_al+tlrc';
epiFilename_2='smnudsbold_e2_tlrc_al+tlrc';
epiFilename_3='smnudsbold_e3_tlrc_al+tlrc';
biopacFilename='biopac.mat';
TEs=[12.8 34.13 55.46];
nEchos=numel(TEs);
TR=2.8;
alpha=0.05;


allOnsets=zeros(nSubjects,2);
allBoldsRoi=cell(nSubjects,1);
allBoldsControlRoi=cell(nSubjects,1);


%%
for s=1:1%nSubjects
    s
    path=['../data/' subjStr{s} '/NII/'];
    biopacPath=['../data/' subjStr{s} '/BIOPAC/'];
    
    [err, roiMask, InfoMask, ~] = BrikLoad (fullfile(path,maskFilename));
    if err || useAverageROI
        [~, roiMask, InfoMask, ~] = BrikLoad (fullfile(path,altMaskFilename));
    end
    
    [err, controlRoiMask, InfoMask, ~] = BrikLoad (fullfile(path,controlMaskFilename));
    if err || useAverageROI
        [~, controlRoiMask, InfoMask, ~] = BrikLoad (fullfile(path,altControlMaskFilename));
    end
    
    [~, brainMask, Info, ~] = BrikLoad (fullfile(path,brainMaskFilename));
    [~, bold_1, Info, ~] = BrikLoad (fullfile(path,epiFilename_1));
    [~, bold_2, Info, ~] = BrikLoad (fullfile(path,epiFilename_2));
    [~, bold_3, Info, ~] = BrikLoad (fullfile(path,epiFilename_3));
    
    % compute the mask
    finalMask=roiMask>0 & brainMask>0; % in the laser beam AND in the brain mask
    finalControlMask=controlRoiMask>0 & brainMask>0;
    
    maskInds=find(finalMask);
    controlMaskInds=find(finalControlMask);
    
    % grab the bold in the ROI
    tmp=permute(bold_1,[4 1 2 3]);
    bold(1,:,:)=tmp(:,:);
    
    tmp=permute(bold_2,[4 1 2 3]);
    bold(2,:,:)=tmp(:,:);
    
    tmp=permute(bold_3,[4 1 2 3]);
    bold(3,:,:)=tmp(:,:);
    
    boldMasked=bold(:,:,finalMask);
    allBoldsRoi{s,1}=boldMasked;
    
    boldControlMasked=bold(:,:,finalControlMask);
    allBoldsControlRoi{s,1}=boldControlMasked;
    
    % laser onset times from biopac
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
    
end


%%
% statistical testing
allBolds=cat(3,allBoldsRoi{:});
nRoiVoxels=size(allBolds,3);
arModelOrder=7;
preLaserTR=214;
for e=1:nEchos
    for v=1:nRoiVoxels
%         tb=squeeze(allBolds(e,1:preLaserTR,v));
%         tb=tb(:);
%         tpl=tplitz(tb,arModelOrder);
%         tpl=tpl(:,2:end-1);
%         h=pinv(tpl)*tb;
%         tbf=filter(1,h,[tb(50:50+arModelOrder-1); zeros(nTRs-arModelOrder,1)]);
        
        tb=squeeze(allBolds(e,1:preLaserTR,:));
        varModel = var(tb,arModelOrder);
        testb=sigma*randn(arModelOrder,1);
        arPred= forecast(arModel,testb,nTRs);
    end
end



%%
[bf,af]=butter(4,0.5,'low');
allBolds=filter(bf,af,allBolds,[],2);
allControlBolds=filter(bf,af,allControlBolds,[],2);

% normalize bolds to 0 at start time of laser
muPre=mean(allBolds(:,TRon,:),2);
allBoldsNorm=allBolds-repmat(muPre,[1 nTRs 1]);

muPreControl=mean(allControlBolds(:,TRon,:),2);
allControlBoldsNorm=allControlBolds-repmat(muPreControl,[1 nTRs 1]);

[sems,mus,stds]=nansem(allBoldsNorm,3);
[semsc,musc,stdsc]=nansem(allControlBoldsNorm,3);

%%
%grand mean in roi
nTRs=size(allBolds,2);
pval=NaN(nEchos,nTRs-TRon+1);
alpha=0.01;
for e=1:nEchos
    for t=TRon:nTRs
        pval(e,t)=ranksum(squeeze(allBolds(e,t,:)),squeeze(allControlBolds(e,t,:)),'tail','right');
    end
end
[pfdr,isSig]=fdr(pval,alpha);

% figure
time=(0:644)*TR;
tLaser=[600 1200];

figure;
for e=1:nEchos
    subplot(nEchos,1,e); hold on
    shadedErrorBar(time,mus(e,:),sems(e,:),'b',1); title('Echo 1');
    shadedErrorBar(time,musc(e,:),semsc(e,:),'r',1); title('Echo 1');     yl=ylim;
    plot(time(isSig(e,:)),yl(1)*ones(sum(isSig(e,:)),1),'*k','LineStyle','none');

    harea=area(tLaser,yl(2)*[1 1],'BaseValue',yl(1));
    set(harea,'FaceColor',[0.7 0.7 0.7]);
    set(harea,'FaceAlpha',0.25);
    set(harea,'EdgeColor','none');
end

%%
% grand mean in control roi

figure;
subplot(411); hold on
plot(time,grandControlBold(1,:)); title('Echo 1');
yl=ylim;
harea=area(tLaser,yl(2)*[1 1],'BaseValue',yl(1));
set(harea,'FaceColor',[0.7 0.7 0.7]);
set(harea,'FaceAlpha',0.5);
set(harea,'EdgeColor','none');
subplot(412); hold on
plot(time,grandControlBold(2,:)); title('Echo 2');
yl=ylim;
harea=area(tLaser,yl(2)*[1 1],'BaseValue',yl(1));
set(harea,'FaceColor',[0.7 0.7 0.7]);
set(harea,'FaceAlpha',0.5);
set(harea,'EdgeColor','none');
subplot(413); hold on
plot(time,grandControlBold(3,:)); title('Echo 3');
yl=ylim;
harea=area(tLaser,yl(2)*[1 1],'BaseValue',yl(1));
set(harea,'FaceColor',[0.7 0.7 0.7]);
set(harea,'FaceAlpha',0.5);
set(harea,'EdgeColor','none');
subplot(414); hold on
plot(time,mean(grandControlBold)); title('Echo mean');
yl=ylim;
harea=area(tLaser,yl(2)*[1 1],'BaseValue',yl(1));
set(harea,'FaceColor',[0.7 0.7 0.7]);
set(harea,'FaceAlpha',0.5);
set(harea,'EdgeColor','none');
print('-dpng','-r600',['../figures/grandMeanControltimeCourse ' date]);

%%

%% 
% % show subject means
% sBold=zeros(nEchos,nTRs,nSubjects);
% for s=1:nSubjects, sBold(:,:,s)=mean(allBoldsRoi{s},3); end
% musBold=mean(sBold,3);
% figure
% for e=1:nEchos
%     subplot(nEchos,1,e);
%     plot(squeeze(sBold(e,:,:)));
%     hold on
%     plot(musBold(e,:),'k','LineWidth',4);
% end
% 
% %%
% figure
% for s=1:nSubjects
%     plot(sBold(2,:,s),'k');
%     title(subjStr{s});
%     pause
% end
