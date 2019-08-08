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
hs(1)=subplot(3,1,1); hold on
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
muBetasPos=cellfun(@mean,pscposb)*100;
sigmaBetasPos=cellfun(@std,pscposb)*100;
muBetasNeg=cellfun(@mean,pscnegb)*100;
sigmaBetasNeg=cellfun(@std,pscnegb)*100;
hs(2)=subplot(3,1,2); hold on
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

set(hs(2),'YLim',[-20 20]);
xl=xlim; yl=ylim;
harea1=area([xl(1) xl(2)],[yl(2) yl(2)],'BaseValue',0,'FaceColor',[0.6 0.6 0.6],'FaceAlpha',0.15,'LineStyle','none');
harea2=area([xl(1) xl(2)],[yl(1) yl(1)],'BaseValue',0,'FaceColor',[0.3 0.3 0.3],'FaceAlpha',0.15,'LineStyle','none');

xl=get(hs(1),'xlim'); set(hs(2),'XLim',xl);

set(get(hs(2),'Ylabel'),'String','% Change');
%set(get(hs(2),'Ylabel'),'Interpreter','Latex');
set(hs(2),'Xtick',1:20);
set(hs(2),'XTickLabel',{'S1','S2','S3','S4','S5', ...
    'S6','S7','S8','S9','S10', ...
    'S11','S12','S13','S14','S15', ...
    'S16','S17','S18','S19','S20'});
set(hs(2),'XTickLabelRotation',45);
set(hs(2),'Ytick',-20:20:20);

set(hbar(4),'FaceColor',get(hbar(1),'FaceColor'));
set(hbar(5),'FaceColor',get(hbar(2),'FaceColor'));
set(hbar(6),'FaceColor',get(hbar(3),'FaceColor'));

%% TIME COURSE
% %%
% grand average signal
%hs(nSubjShow+2)=subplot(5,2,10); hold on;
hs(5)=subplot(3,2,5); hold on
for e=1:3
    plot(time,grandMeanPos(:,e),'LineWidth',1.5);
end
axis tight
set(get(hs(5),'Xlabel'),'String','Time (mm:ss)');
set(hs(5),'Xtick',[0 600 1200 1800]);
set(hs(5),'XtickLabel',{'0:00','10:00','20:00','30:00'});
yl=ylim;
set(hs(5),'YTick',-0.2:0.2:0.2);
plot([600 600],[yl(1) yl(2)],'-k');
plot([1200 1200],[yl(1) yl(2)],'-k');
harea=area([600 1200],[yl(2) yl(2)],'BaseValue',yl(1),'FaceColor',[0.7 0.7 0.7],'FaceAlpha',0.15,'LineStyle','none');
%htit(5)=title('\beta>0','FontWeight','normal');
set(get(hs(5),'ylabel'),'String','BOLD (a.u.)');
htext3=text(30,0.25,'$\beta>0$');set(htext3,'Interpreter','Latex');set(htext3,'FontSize',14);
htext3=text(650,-0.2,'Light on'); set(htext3,'FontSize',14);


hs(6)=subplot(3,2,6); hold on
for e=1:3
    plot(time,grandMeanNeg(:,e),'LineWidth',1.5);
end
axis tight
set(get(hs(6),'Xlabel'),'String','Time (mm:ss)');
set(hs(6),'Xtick',[0 600 1200 1800]);
set(hs(6),'XtickLabel',{'0:00','10:00','20:00','30:00'});
set(hs(6),'ylim',yl);
set(hs(6),'YTick',-0.2:0.2:0.2);
set(hs(6),'YTickLabel',{'','',''});
plot([600 600],[yl(1) yl(2)],'-k');
plot([1200 1200],[yl(1) yl(2)],'-k');
harea=area([600 1200],[yl(2) yl(2)],'BaseValue',yl(1),'FaceColor',[0.7 0.7 0.7],'FaceAlpha',0.15,'LineStyle','none');
%htit(5)=title('\beta<0','FontWeight','normal');
htext4=text(30,0.25,'$\beta<0$');set(htext4,'Interpreter','Latex');set(htext4,'FontSize',14);
htext4=text(650,0.22,'Light on'); set(htext4,'FontSize',14);

%%
hins = axes('position', [0.14 0.865 0.2 0.05]); hold on
plot(hins,1:999,[zeros(1,333) ones(1,333) zeros(1,333)],'k','LineWidth',2); box off
htxt=text(360,1.37,'Design');
set(hins,'XTick',[]); set(hins,'YTick',[]);
xl=xlim; yl=ylim;
harea3=area([xl(1) xl(2)],[yl(2) yl(2)],'BaseValue',0,'FaceColor',[0.6 0.6 0.6],'FaceAlpha',0.15,'LineStyle','none');

sublabel([hs(1) hs(2) hs(5) hs(6)],0,-45,'FontWeight','Bold','FontSize',16);
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

