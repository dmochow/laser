% ROI-based analysis
% look at percentage of significant voxels, beta weights, time course
% at both true and control rois
%
% examine scaling of all BOLD series (output of 3dTproject)
clear all; close all; clc
addpath(genpath('~/PROJECTS/COMMON'));

useAverageROI=0; % if we want to apply the same average ROI to all subjects
subjStr={'S02','S03','S04','S05','S06','S07','S09','S10','S11','S12','S13','S14','S15'};
nSubjects=numel(subjStr);
writeBriks=1;

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

%allMaskInds=cell(nSubjects,1);
allOnsets=zeros(nSubjects,2);
isSigAll=cell(nSubjects,1);
isSigAllControl=cell(nSubjects,1);
allBetas=cell(nSubjects,1);
allBetasControl=cell(nSubjects,1);
allBoldsRoi=cell(nSubjects,1);
allBoldsControlRoi=cell(nSubjects,1);
nSig=zeros(nEchos,nSubjects);
propSig=zeros(nEchos,nSubjects);
nSigControl=zeros(nEchos,nSubjects);
propSigControl=zeros(nEchos,nSubjects);

%%
for s=1:nSubjects
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

    
    % start GLM here (real ROI)
    % for each voxel in ROI, convert to 1D time series
    nTRs=size(boldMasked,2);
    nROIvoxels=size(boldMasked,3);
    nControlROIvoxels=size(boldControlMasked,3);
    pval=zeros(nEchos,nROIvoxels);
    pvalControl=zeros(nEchos,nControlROIvoxels);
    betas=zeros(nEchos,nROIvoxels,nTRsLaser);
    betasControl=zeros(nEchos,nControlROIvoxels,nTRsLaser);

    for e=1:nEchos
        % true ROI
        for v=1:nROIvoxels
            Y=squeeze(boldMasked(e,:,v));
            Y=Y(:);
            
            % FIR
            Xf=eye(nTRs);
            Xf=Xf(:,onsetTimeTR:offsetTimeTR);
            
            % boxcar
%             Xf=zeros(nTRs,1);
%             Xf(onsetTimeTR:offsetTimeTR)=1;
            
            statsf=myGLM(Y,Xf);
            ssef=statsf.sse;
            doff=statsf.dof;
            betas(e,v,:)=statsf.B(1:end-1);
            
            % reduced model
            %Xr=ones(nTRs,1);
            Xr=[]; % myGLM adds in ones column
            statsr=myGLM(Y,Xr);
            sser=statsr.sse;
            dofr=statsr.dof;
            Fstat=((sser-ssef)/(dofr-doff)) / (ssef/doff);
            pval(e,v)=fcdf(Fstat,dofr-doff,doff,'upper');
        end
        
        % control ROI
        for v=1:nControlROIvoxels
            Yc=squeeze(boldControlMasked(e,:,v));
            Yc=Yc(:);
            statsfc=myGLM(Yc,Xf);
            ssefc=statsfc.sse;
            doffc=statsfc.dof;
            betasControl(e,v,:)=statsfc.B(1:end-1);
            
            % reduced model
            statsrc=myGLM(Yc,Xr);
            sserc=statsrc.sse;
            dofrc=statsrc.dof;
            Fstatc=((sserc-ssefc)/(dofrc-doffc)) / (ssefc/doffc);
            pvalControl(e,v)=fcdf(Fstatc,dofrc-doffc,doffc,'upper');
        end
    end
    allBetas{s}=betas;
    allBetasControl{s}=betasControl;
    
    % correct for multiple comparisons
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
    
    if writeBriks
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
    
    [p_fdr_control, p_masked_control] = fdr( pvalControl, alpha);
    isSig_fdr_control=p_masked_control;
    
    isSig_e1_control=zeros(size(finalControlMask));
    isSig_e1_control(controlMaskInds)=isSig_fdr_control(1,:);
    
    isSig_e2_control=zeros(size(finalControlMask));
    isSig_e2_control(controlMaskInds)=isSig_fdr_control(2,:);
    
    isSig_e3_control=zeros(size(finalControlMask));
    isSig_e3_control(controlMaskInds)=isSig_fdr_control(3,:);
    
    isSigAllControl{s,1}=isSig_e1_control;
    isSigAllControl{s,2}=isSig_e2_control;
    isSigAllControl{s,3}=isSig_e3_control;
    
    nSig(:,s)=sum(p_masked,2);
    nSigControl(:,s)=sum(p_masked_control,2);
    propSig(:,s)=mean(p_masked,2);
    propSigControl(:,s)=mean(p_masked_control,2);
    
end

%%
% make bar graph from proportion of significant voxels
xvals=[1 2 3];
deloff=0.15;
[sems,mus] = nansem( propSig,2 );
[semsc,musc] = nansem( propSigControl,2 );
figure; hold on
hbar=bar(xvals,[mus musc]);
set(hbar(1),'EdgeColor','none');
set(hbar(1),'BarWidth',0.9);
set(hbar(1),'FaceAlpha',0.9);
set(hbar(2),'EdgeColor','none');
set(hbar(2),'FaceAlpha',0.9);
herb=errorbar(xvals-deloff,mus,sems,'LineStyle','none');
herbc=errorbar(xvals+deloff,musc,semsc,'LineStyle','none');
set(herb,'Color','k');
set(herbc,'Color','k');
xlabel('Echo')
ylabel('Percent Significant Voxels')
hlg=legend('ROI','contralateral');
set(hlg,'box','off');
set(hlg,'Location','northwest');
set(gca,'XTick',xvals);
set(gca,'YTick',0:0.05:0.25);
set(gca,'YTickLabel',{'0','5','10','15','20','25'});

print('-dpng',['../figures/percentSigVoxels-' date]);

%%
% grand mean in roi
subs=setdiff(1:nSubjects,[]);
time=(0:644)*TR;
tLaser=[600 1200];
grandBold=mean(cat(3,allBoldsRoi{subs}),3);
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
print('-dpng','-r600',['../figures/grandMeanROItimeCourse ' date]);

%%
% because I'm eager
time=(0:644)*TR;
tLaser=[600 1200];
grandControlBold=mean(cat(3,allBoldsControlRoi{:}),3);
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
