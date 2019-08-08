% 05.22.19
% separated positive from negative activations
% computed psc

% 06.24.19
% adding ROI drawing + design
% removing time courses

clear all; close all; clc
PRINT_FIGURE=1;

savedDataFilename='../data/precomputed/revFig1data.mat';
subjStrs={'S02','S03','S04','S05','S06','S07','S09','S10','S11','S12','S13'...
    ,'S14','S15','S16','S17','S18','S19','S20','S21','S22'};
nSubjects=numel(subjStrs);
TR=2.8;
nTRs=645; % hard-coded but needs to be changed if this script is to be adapted for the toggling study
time=(0:nTRs-1)*TR;
alpha=0.05; % significance level
gmIndex=2; % this will never change
roiMaskFilename='resampled_roi_r21_z39+tlrc.BRIK';
altRoiMaskFilename='resampled_mu_roi_r21_z39+tlrc.BRIK'; % no markers
segFilename='resampled_Classes+tlrc.BRIK';
boldPrefix='smoothedOutputBoldEcho'; % use this to test different preprocs
biopacFilename='biopac.mat';
nSubjShow=7;
segs=[1:215;216:430;431:645];


load(savedDataFilename,'r','r2','allts','issig','allBetas','allOnsets');

%%
% construct grand mean time course for significant voxels
% split based on positive or negative beta
muts=zeros(nTRs,3);
allsigposts=cell(3,1);
allsignegts=cell(3,1);
for s=1:nSubjects
    for e=1:3
        % shift this subject's data to enforce laser on at TR=215
        onsetTR=floor(allOnsets(s,1)/TR);
        tshift=215-onsetTR;
        allts{e,s}=circshift(allts{e,s},[tshift 0]);
        allsigposts{e}=cat(2,allsigposts{e},allts{e,s}(:,issig{e,s} & allBetas{e,s}(1,:)>0));
        allsignegts{e}=cat(2,allsignegts{e},allts{e,s}(:,issig{e,s} & allBetas{e,s}(1,:)<0));
        
        sigposts{e,s}=allts{e,s}(:,issig{e,s} & allBetas{e,s}(1,:)>0);
        signegts{e,s}=allts{e,s}(:,issig{e,s} & allBetas{e,s}(1,:)<0);
        muposts(:,e,s)=mean(sigposts{e,s},2);
        munegts(:,e,s)=mean(signegts{e,s},2);
        
        % figure out positive or negative betas
        % note: using B(1) here (acute effect)
        % B(2) is after effect
        % B(3) is constant term
        tbetas=allBetas{e,s}(1,:);
        tbetas(~issig{e,s})=0; % remove non-significant ones
        propSigPos(e,s)= mean(tbetas>0);
        propSigNeg(e,s)= mean(tbetas<0);
        sumSigPos(e,s)= numel(find(tbetas>0));
        sumSigNeg(e,s)= numel(find(tbetas<0));
        
        sigbetaspos3=allBetas{e,s}(:,find(issig{e,s} & allBetas{e,s}(1,:)>0   ));
        pscpos{e,s}=100*abs(sigbetaspos3(1,:)./sigbetaspos3(3,:));
        pscposb{e,s}=sigbetaspos3(1,:);
        %pscposb{e,s}=100*sigbetaspos3(1,:)./ abs(sigbetaspos3(3,:));
        
        sigbetasneg3=allBetas{e,s}(:,find(issig{e,s} & allBetas{e,s}(1,:)<0   ));
        pscneg{e,s}=100*abs(sigbetasneg3(1,:)./sigbetasneg3(3,:));
        pscnegb{e,s}=sigbetasneg3(1,:);
        %pscnegb{e,s}=100*sigbetasneg3(1,:)./abs(sigbetasneg3(3,:));
        
    end
end

%%
% mean ts of each subject
grandMeanPos=nanmean(muposts,3);
grandMeanNeg=nanmean(munegts,3);
grandSigmaPos=nanstd(muposts,[],3);
grandSigmaNeg=nanstd(munegts,[],3);
tb=round(600/TR); ta=round(800/TR);
[grandMeanPos(tb) grandSigmaPos(tb) grandMeanPos(ta) grandSigmaPos(ta)]
0.5*(grandSigmaPos(tb)+grandSigmaPos(ta))

%% # SIGNIFICANT
figure
hs(1)=subplot(3,1,2); hold on
hbar=bar(1:20,[100*propSigPos' -100*propSigNeg']);
xl=xlim; yl=ylim;
harea1=area([xl(1) xl(2)],[yl(2) yl(2)],'BaseValue',0,'FaceColor',[0.6 0.6 0.6],'FaceAlpha',0.15,'LineStyle','none');
harea2=area([xl(1) xl(2)],[yl(1) yl(1)],'BaseValue',0,'FaceColor',[0.3 0.3 0.3],'FaceAlpha',0.15,'LineStyle','none');
htext1=text(1,30,'$\beta>0$');set(htext1,'Interpreter','Latex');set(htext1,'FontSize',14);
htext2=text(1,-30,'$\beta<0$');set(htext2,'Interpreter','Latex');set(htext2,'FontSize',14);
set(get(hs(1),'ylabel'),'String','% sig voxels');
set(hs(1),'Xtick',1:20);
set(hs(1),'XTickLabel',{'S1','S2','S3','S4','S5', ...
    'S6','S7','S8','S9','S10', ...
    'S11','S12','S13','S14','S15', ...
    'S16','S17','S18','S19','S20'});
set(hs(1),'XTickLabelRotation',45);
set(hs(1),'Ytick',[-50 0 50 100]);
set(hs(1),'YTickLabel',{'50','0','50','100'});
hlg0=legend('Echo 1','Echo 2','Echo 3');
set(hlg0,'Orientation','horizontal');
set(hlg0,'box','off');
set(hbar(4),'FaceColor',get(hbar(1),'FaceColor'));
set(hbar(5),'FaceColor',get(hbar(2),'FaceColor'));
set(hbar(6),'FaceColor',get(hbar(3),'FaceColor'));

%% EFFECT SIZE
muBetasPos=cellfun(@mean,pscposb)*1;
sigmaBetasPos=cellfun(@std,pscposb)*1;
muBetasNeg=cellfun(@mean,pscnegb)*1;
sigmaBetasNeg=cellfun(@std,pscnegb)*1;

hs(2)=subplot(3,1,3); hold on
hbar=bar(1:20,[muBetasPos' muBetasNeg']);

ngroups=20;
nbars=6;
groupwidth = min(0.8, nbars/(nbars + 1.5));
y=[muBetasPos' muBetasNeg'];
err=[sigmaBetasPos' sigmaBetasNeg'];
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    herr(:,i)=errorbar(x,y(:,i), err(:,i), '.k','CapSize',0);
end
%hold off

%set(hs(2),'YLim',[-20 20]);
xl=xlim; yl=ylim;
harea1=area([xl(1) xl(2)],[yl(2) yl(2)],'BaseValue',0,'FaceColor',[0.6 0.6 0.6],'FaceAlpha',0.15,'LineStyle','none');
harea2=area([xl(1) xl(2)],[yl(1) yl(1)],'BaseValue',0,'FaceColor',[0.3 0.3 0.3],'FaceAlpha',0.15,'LineStyle','none');

xl=get(hs(1),'xlim'); set(hs(2),'XLim',xl);

%set(get(hs(2),'Ylabel'),'String','% Change');
set(get(hs(2),'Ylabel'),'String','Beta weights');
%set(get(hs(2),'Ylabel'),'Interpreter','Latex');
set(hs(2),'Xtick',1:20);
set(hs(2),'XTickLabel',{'S1','S2','S3','S4','S5', ...
    'S6','S7','S8','S9','S10', ...
    'S11','S12','S13','S14','S15', ...
    'S16','S17','S18','S19','S20'});
set(hs(2),'XTickLabelRotation',45);
%set(hs(2),'Ytick',-20:20:20);

set(hbar(4),'FaceColor',get(hbar(1),'FaceColor'));
set(hbar(5),'FaceColor',get(hbar(2),'FaceColor'));
set(hbar(6),'FaceColor',get(hbar(3),'FaceColor'));

% %% TIME COURSE
% % %%
% % grand average signal
% %hs(nSubjShow+2)=subplot(5,2,10); hold on;
% hs(5)=subplot(3,2,5); hold on
% for e=1:3
%     plot(time,grandMeanPos(:,e),'LineWidth',1.5);
% end
% axis tight
% set(get(hs(5),'Xlabel'),'String','Time (mm:ss)');
% set(hs(5),'Xtick',[0 600 1200 1800]);
% set(hs(5),'XtickLabel',{'0:00','10:00','20:00','30:00'});
% yl=ylim;
% set(hs(5),'YTick',-0.2:0.2:0.2);
% plot([600 600],[yl(1) yl(2)],'-k');
% plot([1200 1200],[yl(1) yl(2)],'-k');
% harea=area([600 1200],[yl(2) yl(2)],'BaseValue',yl(1),'FaceColor',[0.7 0.7 0.7],'FaceAlpha',0.15,'LineStyle','none');
% %htit(5)=title('\beta>0','FontWeight','normal');
% set(get(hs(5),'ylabel'),'String','BOLD (a.u.)');
% htext3=text(30,0.25,'$\beta>0$');set(htext3,'Interpreter','Latex');set(htext3,'FontSize',14);
% htext3=text(650,-0.2,'Light on'); set(htext3,'FontSize',14);
% 
% 
% hs(6)=subplot(3,2,6); hold on
% for e=1:3
%     plot(time,grandMeanNeg(:,e),'LineWidth',1.5);
% end
% axis tight
% set(get(hs(6),'Xlabel'),'String','Time (mm:ss)');
% set(hs(6),'Xtick',[0 600 1200 1800]);
% set(hs(6),'XtickLabel',{'0:00','10:00','20:00','30:00'});
% set(hs(6),'ylim',yl);
% set(hs(6),'YTick',-0.2:0.2:0.2);
% set(hs(6),'YTickLabel',{'','',''});
% plot([600 600],[yl(1) yl(2)],'-k');
% plot([1200 1200],[yl(1) yl(2)],'-k');
% harea=area([600 1200],[yl(2) yl(2)],'BaseValue',yl(1),'FaceColor',[0.7 0.7 0.7],'FaceAlpha',0.15,'LineStyle','none');
% %htit(5)=title('\beta<0','FontWeight','normal');
% htext4=text(30,0.25,'$\beta<0$');set(htext4,'Interpreter','Latex');set(htext4,'FontSize',14);
% htext4=text(650,0.22,'Light on'); set(htext4,'FontSize',14);

%%
nShow=3;
subj2show=[7;13;19];
zcoors=[98;93;107];
%zcoor=110; % S20

for s=1:nShow

    zcoor=zcoors(s);
    pathToData=['../data/myoutput/' subjStrs{subj2show(s)} '/'];
    [~, anat, ~, ~] = BrikLoad (fullfile(pathToData,'anat+tlrc'));
    anat=anat(end:-1:1,:,:);
    
    anatRoiMaskFilename='roi_r21_z39+tlrc.BRIK';
    [~, anatRoiMask, ~, ~] = BrikLoad (fullfile(pathToData,anatRoiMaskFilename));
    anatRoiMask=anatRoiMask(end:-1:1,:,:);
    anatRoiMask=logical(anatRoiMask);
    
    anatBrainMaskFilename='brain_mask+tlrc.BRIK';
    [~, anatBrainMask, ~, ~] = BrikLoad (fullfile(pathToData,anatBrainMaskFilename));
    anatBrainMask=anatBrainMask(end:-1:1,:,:);
    
    issigMask=zeros(size(anatBrainMask));
    
    colorMask=anatRoiMask(:,:,zcoor)>0 & anatBrainMask(:,:,zcoor)>0;
    
    hs(3+s)=subplot(3,6,s);
    %h = colorBrain(anat(:,:,zcoor),anatRoiMask(:,:,zcoor),issigMask(:,:,zcoor)==1,anatBrainMask(:,:,zcoor),1);
    h = colorBrain(anat(:,:,zcoor),colorMask,issigMask(:,:,zcoor)==1,anatBrainMask(:,:,zcoor),1);
    freezeColors
    set(h,'FaceAlpha',1);
    htit(s)=title(['Subject ' num2str(subj2show(s))],'FontWeight','normal')
    moveTitle(htit(s),0,210,0);
end

%%
%hins = axes('position', [0.14 0.865 0.2 0.05]); hold on
hs(3)=subplot(322);  hold on
%plot(1:999,[zeros(1,333) ones(1,333) zeros(1,333)],'k','LineWidth',2); box off
%plot(1:999,[zeros(1,333) zeros(1,333) ones(1,333)],'r','LineWidth',2); 
designMatrix(:,1)=[zeros(1,333) ones(1,333) zeros(1,333)];
designMatrix(:,2)=[zeros(1,333) zeros(1,333) ones(1,333)];
himg=imagesc([1 2],[0 1800],designMatrix); 
cm=colormap; colormap(cm([1 end],:));
axis ij;
%htxt=text(360,1.37,'Design');
htits=title('Design Matrix','FontWeight','normal'); 
set(hs(3),'Ytick',[0 600 1200 1800]); 
set(hs(3),'YTickLabel',{'0','10','20','30'}); 
set(hs(3),'XTick',[1.25 1.75]);set(hs(3),'XTickLabel',{'1','2'}); 
%set(hs(3),'XTickLabel',{'0','10','20','30'});
xl=xlim; yl=ylim;
xlim([1 2]); ylim([0 1800]);
ylabel('Time (min)');
xlabel('Column');
subpos=get(hs(3),'Position');
set(hs(3),'Position',[subpos(1)+0.1 subpos(2) subpos(3)*0.5 subpos(4)]);
hcb=colorbar; set(hcb,'Ticks',[0.25 0.75]); set(hcb,'TickLabels',{'0','1'});
%harea3=area([xl(1) xl(2)],[yl(2) yl(2)],'BaseValue',0,'FaceColor',[0.6 0.6 0.6],'FaceAlpha',0.15,'LineStyle','none');

sublabel([hs(4) hs(3) hs(1) hs(2)],-10,-35,'FontWeight','Bold','FontSize',16);
%
if PRINT_FIGURE
    %%
    %
    %  PREFORMATTED
    %  TEXT
    %
    figStr=['../figures/revfig1-' date '.png'];
    print('-dpng',figStr); % ../figures/revfig1
    crop(figStr,0);
end

%%
% stats that are reported in the papar
%
% number of voxels in ROI
%
%
[~,nVoxelsROI]=cellfun(@size,allts);
nVoxelsROI=nVoxelsROI(1,:);
munVoxelsROI=mean(nVoxelsROI)
stdnVoxelsROI=std(nVoxelsROI)
%
% number of subjects with sig voxels
nSigSubjects=numel(find((sum(r2>0,1))))
nSigPosSubjects=numel(find((sum(propSigPos>0,1))))
nSigNegSubjects=numel(find((sum(propSigNeg>0,1))))
nSigBothSubjects=numel(find((sum(propSigNeg>0 & propSigPos>0,1))))

mupctSig=mean(r(:));
sigmapctSig=std(r(:));

%%
muPSC=nanmean(muBetasPos,2)
sigmaPSC=nanstd(muBetasPos,[],2)