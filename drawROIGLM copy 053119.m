%% 05.22.19

% separated positive from negative activations
% computed psc

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
        
        sigbetasneg3=allBetas{e,s}(:,find(issig{e,s} & allBetas{e,s}(1,:)<0   ));
        pscneg{e,s}=100*abs(sigbetasneg3(1,:)./sigbetasneg3(3,:));
        pscnegb{e,s}=sigbetasneg3(1,:);
        
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
hs(1)=subplot(3,1,1);
bar(1:20,[100*propSigPos' -100*propSigNeg']);
set(get(hs(1),'ylabel'),'String','% sig voxels');
set(hs(1),'Xtick',1:20);
set(hs(1),'XTickLabel',{'S1','S2','S3','S4','S5', ...
    'S6','S7','S8','S9','S10', ...
    'S11','S12','S13','S14','S15', ...
    'S16','S17','S18','S19','S20'});
set(hs(1),'XTickLabelRotation',45);
hlg0=legend('Echo 1 +','Echo 2 +','Echo 3 +','Echo 1 -','Echo 2 -','Echo 3 -');
    
%% EFFECT SIZE
hs(2)=subplot(3,2,3);
violin(pscposb(1,propSigPos(1,:)>0));
hs(3)=subplot(3,2,4);
violin(pscnegb(1,propSigNeg(1,:)>0));
%ylim([0 300]);


%% TIME COURSE
% %%
% grand average signal
%hs(nSubjShow+2)=subplot(5,2,10); hold on;
hs(5)=subplot(3,2,5); hold on
for e=1:3
    plot(time,grandMeanPos(:,e),'LineWidth',1.5);
    %plot(time,mean(allsigposts{e},2),'LineWidth',1.5);
    %plot(time,mean(allsignegts{e},2),'LineWidth',1.5);
    axis tight
    yl=ylim;
    plot([600 600],[yl(1) yl(2)],'-k');
    plot([1200 1200],[yl(1) yl(2)],'-k');
    harea=area([600 1200],[yl(2) yl(2)],'BaseValue',yl(1),'FaceColor',[0.7 0.7 0.7],'FaceAlpha',0.15,'LineStyle','none');
end
htit(5)=title('Beta Positive');

hs(3)=subplot(3,2,6); hold on
for e=1:3
    plot(time,grandMeanNeg(:,e),'LineWidth',1.5);
    %plot(time,mean(allsigposts{e},2),'LineWidth',1.5);
    %plot(time,mean(allsignegts{e},2),'LineWidth',1.5);
    axis tight
    yl=ylim;
    plot([600 600],[yl(1) yl(2)],'-k');
    plot([1200 1200],[yl(1) yl(2)],'-k');
    harea=area([600 1200],[yl(2) yl(2)],'BaseValue',yl(1),'FaceColor',[0.7 0.7 0.7],'FaceAlpha',0.15,'LineStyle','none');
end
htit(5)=title('Beta Negative');

%% signal powers
pb=mean(grandMeanPos(segs(1,:),:).^2,1);
pd=mean(grandMeanPos(segs(2,:),:).^2,1);
pa=mean(grandMeanPos(segs(3,:),:).^2,1);
ppc=100*(pd-pb)./pb;
ppc2=100*(pa-pb)./pb;

pb=mean(grandMeanNeg(segs(1,:),:).^2,1);
pd=mean(grandMeanNeg(segs(2,:),:).^2,1);
pa=mean(grandMeanNeg(segs(3,:),:).^2,1);
ppc3=100*(pd-pb)./pb;
ppc4=100*(pa-pb)./pb;



% 
% %%
% %
% %figure
% % find best subjects
% %[~,sortind]=sort(mean(r),'descend');
% [~,sortind]=sort(mean(r2),'descend');
% 
% for s=1:nSubjShow
%     %hs(s+1)=subplot(6,2,4+s); hold on
%     hs(s+1)=subplot(5,2,2+s); hold on
%     for e=1:3
%         thisTs=allts{e,sortind(s)};
%         vinds=issig{e,sortind(s)};
%         muThisTs=mean(thisTs(:,vinds==1),2);
%         plot(time,muThisTs);
%     end
%     %ylim([-1 1]);
%     axis tight
%     yl=ylim;
%     plot([allOnsets(s,1) allOnsets(s,1)],[yl(1) yl(2)],'-k');
%     plot([allOnsets(s,2) allOnsets(s,2)],[yl(1) yl(2)],'-k');
%     harea=area([allOnsets(s,1) allOnsets(s,2)],[yl(2) yl(2)],'BaseValue',yl(1),'FaceColor',[0.7 0.7 0.7],'FaceAlpha',0.15,'LineStyle','none');
% 
%     htxt=text(50,yl(1)+0.15,['S' num2str(sortind(s))]);
%     
%     if s==1, 
%         hlg=legend('Echo 1','Echo 2','Echo 3'); 
%         htxt=text(725,yl(1)+0.25,['Light on']);
%     end
%     
% end
% 
% hins = axes('position', [0.14 0.91 0.2 0.05]);
% plot(hins,1:999,[zeros(1,333) ones(1,333) zeros(1,333)],'k','LineWidth',2); box off
% htxt=text(360,1.37,'Design');
% set(hins,'XTick',[]); set(hins,'YTick',[]);
% 
% % grand average signal
% hs(nSubjShow+2)=subplot(5,2,10); hold on;
% for e=1:3
%     plot(time,mean(allsigts{e},2),'LineWidth',1.5);
%     axis tight
%     yl=ylim;
%     plot([600 600],[yl(1) yl(2)],'-k');
%     plot([1200 1200],[yl(1) yl(2)],'-k');
%     harea=area([600 1200],[yl(2) yl(2)],'BaseValue',yl(1),'FaceColor',[0.7 0.7 0.7],'FaceAlpha',0.15,'LineStyle','none');
% end
% htxt=text(50,yl(1)+0.15,'AVG');
% 
% %%
% % esthetics
% set(hs(1),'Xtick',1:20);
% set(hs(1),'Box','off');
% 
% set(hs(2:end),'XLim',[0 1800]);
% set(hs(2:end),'XTick',[0 600 1200 1800]);
% set(hs(2:nSubjShow),'XTickLabel',{'','','',''});
% set(hs(nSubjShow+1),'XTickLabel',{'0:00','10:00','20:00','30:00'});
% set(hs(nSubjShow+2),'XTickLabel',{'0:00','10:00','20:00','30:00'});
% 
% set(hs(2:end), 'YLimSpec', 'Tight');
% 
% set(hs(2:end),'Ytick',[-0.25 0 0.25]);
% set(hs(3:2:end),'YtickLabel',{'','',''});
% %set(hs(3:2:end),'YTickLabel',{'','',''});
% 
% set(get(hs(2),'YLabel'),'String','BOLD (a.u.)');
% set(get(hs(nSubjShow+1),'XLabel'),'String','Time (mm:ss)');
% set(get(hs(nSubjShow+2),'XLabel'),'String','Time (mm:ss)');
% 
% set(hlg0,'box','off','Orientation','horizontal');
% 
% pos=get(hs(1),'Position');
% delx=0; dely=0.05;
% set(hs(1),'Position',[pos(1) pos(2)+dely pos(3) pos(4)]);
% 
% delx=0; dely=-0.025;
% for s=2:nSubjShow+2
%     pos=get(hs(s),'Position');
%     set(hs(s),'Position',[pos(1) pos(2)+dely pos(3) pos(4)]);
% end
% 
% lgPos=get(hlg,'Position');
% delx=-0.09; dely=-0.09;
% set(hlg,'Position',[lgPos(1)+delx lgPos(2)+dely lgPos(3) lgPos(4)]);
% set(hlg,'box','off','Orientation','horizontal');
% 
% sublabel([hs(1) hs(2)],0,-60,'FontWeight','Bold','FontSize',16);
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

% %% 
% % make spatial figure here
%%
% s=20;
% b1=allBetas{1,s}(issig{1,s});
% b2=allBetas{2,s}(issig{1,s});
% b3=allBetas{3,s}(issig{1,s});
% figure;
% hold on
% plot([1 2 3],[b1' b2' b3'],'--o');

% 
% %%
% pathToData=['../data/myoutput/' subjStrs{s} '/'];
% pathToBiopac=['../data/scanner_data/' subjStrs{s} '/'];
% [~, roiMask, info, ~] = BrikLoad (fullfile(pathToData,roiMaskFilename));
% roiMask=logical(roiMask);
% resampledSegFilename=fullfile(pathToData,'resampled_Classes+tlrc.BRIK');
% [~, seg, ~, ~] = BrikLoad (resampledSegFilename);
% greyMask=logical(seg==gmIndex);
% finalMask=roiMask&greyMask;
% betaVol=ts2vol(allBetas{1,s},finalMask);
% 
% %%
% 
% %h = colorBrain(anatShow,isSig1,isSig2,brainMask,crop,nColors)