%% 05.19.19
%
% draw the results of the ROI GLM analysis, forming what should be figure 1
% in the revised manuscript

% move up bottom 3 rows variable amounts

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



load(savedDataFilename,'r','r2','allts','issig','allBetas','allOnsets');

%%
% construct grand mean time course for significant voxels
muts=zeros(nTRs,3);
allsigts=cell(3,1);
for s=1:nSubjects
    for e=1:3
        % shift this subject's data to enforce laser on at TR=215
        onsetTR=floor(allOnsets(s,1)/TR);
        tshift=215-onsetTR;
        allts{e,s}=circshift(allts{e,s},[tshift 0]);
        allsigts{e}=cat(2,allsigts{e},allts{e,s}(:,issig{e,s}));
    end
end

%%
figure
hs(1)=subplot(411);
bar(1:20,100*r');
set(get(hs(1),'ylabel'),'String','% sig voxels');
%set(get(hs(1),'xlabel'),'String','Subject');
set(hs(1),'Xtick',1:20);
set(hs(1),'XTickLabel',{'S1','S2','S3','S4','S5', ...
    'S6','S7','S8','S9','S10', ...
    'S11','S12','S13','S14','S15', ...
    'S16','S17','S18','S19','S20'});
set(hs(1),'XTickLabelRotation',45);
hlg0=legend('Echo 1','Echo 2','Echo 3'); 

%%
% grand average signal
%hs(nSubjShow+2)=subplot(5,2,10); hold on;
hs(2)=subplot(4,1,2);
for e=1:3
    plot(time,mean(allsigts{e},2),'LineWidth',1.5);
    axis tight
    yl=ylim;
    plot([600 600],[yl(1) yl(2)],'-k');
    plot([1200 1200],[yl(1) yl(2)],'-k');
    harea=area([600 1200],[yl(2) yl(2)],'BaseValue',yl(1),'FaceColor',[0.7 0.7 0.7],'FaceAlpha',0.15,'LineStyle','none');
end
htxt=text(50,yl(1)+0.15,'AVG');

%%
%
%figure
% find best subjects
%[~,sortind]=sort(mean(r),'descend');
[~,sortind]=sort(mean(r2),'descend');

for s=1:nSubjShow
    %hs(s+1)=subplot(6,2,4+s); hold on
    hs(s+1)=subplot(5,2,2+s); hold on
    for e=1:3
        thisTs=allts{e,sortind(s)};
        vinds=issig{e,sortind(s)};
        muThisTs=mean(thisTs(:,vinds==1),2);
        plot(time,muThisTs);
    end
    %ylim([-1 1]);
    axis tight
    yl=ylim;
    plot([allOnsets(s,1) allOnsets(s,1)],[yl(1) yl(2)],'-k');
    plot([allOnsets(s,2) allOnsets(s,2)],[yl(1) yl(2)],'-k');
    harea=area([allOnsets(s,1) allOnsets(s,2)],[yl(2) yl(2)],'BaseValue',yl(1),'FaceColor',[0.7 0.7 0.7],'FaceAlpha',0.15,'LineStyle','none');

    htxt=text(50,yl(1)+0.15,['S' num2str(sortind(s))]);
    
    if s==1, 
        hlg=legend('Echo 1','Echo 2','Echo 3'); 
        htxt=text(725,yl(1)+0.25,['Light on']);
    end
    
end

hins = axes('position', [0.14 0.91 0.2 0.05]);
plot(hins,1:999,[zeros(1,333) ones(1,333) zeros(1,333)],'k','LineWidth',2); box off
htxt=text(360,1.37,'Design');
set(hins,'XTick',[]); set(hins,'YTick',[]);

% grand average signal
hs(nSubjShow+2)=subplot(5,2,10); hold on;
for e=1:3
    plot(time,mean(allsigts{e},2),'LineWidth',1.5);
    axis tight
    yl=ylim;
    plot([600 600],[yl(1) yl(2)],'-k');
    plot([1200 1200],[yl(1) yl(2)],'-k');
    harea=area([600 1200],[yl(2) yl(2)],'BaseValue',yl(1),'FaceColor',[0.7 0.7 0.7],'FaceAlpha',0.15,'LineStyle','none');
end
htxt=text(50,yl(1)+0.15,'AVG');

%%
% esthetics
set(hs(1),'Xtick',1:20);
set(hs(1),'Box','off');

set(hs(2:end),'XLim',[0 1800]);
set(hs(2:end),'XTick',[0 600 1200 1800]);
set(hs(2:nSubjShow),'XTickLabel',{'','','',''});
set(hs(nSubjShow+1),'XTickLabel',{'0:00','10:00','20:00','30:00'});
set(hs(nSubjShow+2),'XTickLabel',{'0:00','10:00','20:00','30:00'});

set(hs(2:end), 'YLimSpec', 'Tight');

set(hs(2:end),'Ytick',[-0.25 0 0.25]);
set(hs(3:2:end),'YtickLabel',{'','',''});
%set(hs(3:2:end),'YTickLabel',{'','',''});

set(get(hs(2),'YLabel'),'String','BOLD (a.u.)');
set(get(hs(nSubjShow+1),'XLabel'),'String','Time (mm:ss)');
set(get(hs(nSubjShow+2),'XLabel'),'String','Time (mm:ss)');

set(hlg0,'box','off','Orientation','horizontal');

pos=get(hs(1),'Position');
delx=0; dely=0.05;
set(hs(1),'Position',[pos(1) pos(2)+dely pos(3) pos(4)]);

delx=0; dely=-0.025;
for s=2:nSubjShow+2
    pos=get(hs(s),'Position');
    set(hs(s),'Position',[pos(1) pos(2)+dely pos(3) pos(4)]);
end

lgPos=get(hlg,'Position');
delx=-0.09; dely=-0.09;
set(hlg,'Position',[lgPos(1)+delx lgPos(2)+dely lgPos(3) lgPos(4)]);
set(hlg,'box','off','Orientation','horizontal');

sublabel([hs(1) hs(2)],0,-60,'FontWeight','Bold','FontSize',16);

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