% 07.16.19
% examine state-dependence of BOLD response

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
colors=[0 0.4470 0.7410;0.8500    0.3250    0.0980;0.9290    0.6940    0.1250];

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

[grandSemPos,grandMuPos]=nansem(muposts,3);
[grandSemNeg,grandMuNeg]=nansem(munegts,3);

tb=round(600/TR); ta=round(800/TR);
[grandMeanPos(tb) grandSigmaPos(tb) grandMeanPos(ta) grandSigmaPos(ta)]
0.5*(grandSigmaPos(tb)+grandSigmaPos(ta))

%%
% look at slopes in the baseline period
bpos=NaN(3,nSubjects);
bneg=NaN(3,nSubjects);
mupscpos=NaN(3,nSubjects);
mupscneg=NaN(3,nSubjects);
for e=1:3
    for s=1:nSubjects
        
        y=muposts(segs(1,:),e,s);
        x=[segs(1,:)' ones(size(segs,2),1)];
        if ~sum(isnan(y))
            b=regress(y,x);
            bpos(e,s)=b(1);
            mupscpos(e,s)=mean(pscposb{e,s});
        end
        
        y=munegts(segs(1,:),e,s);
        x=[segs(1,:)' ones(size(segs,2),1)];
        if ~sum(isnan(y))
            b=regress(y,x);
            bneg(e,s)=b(1);
            mupscneg(e,s)=mean(pscnegb{e,s});
        end

    end
end

[sembpos,mubpos]=nansem(bpos,2);
[sembneg,mubneg]=nansem(bneg,2);

% test for significant difference from zero
for e=1:3
    [hh,pvalpos(e)]=ttest(bpos(e,~isnan(bpos(e,:))));
    [hh,pvalneg(e)]=ttest(bneg(e,~isnan(bneg(e,:))));
end

%% TIME COURSE
nRows=3;
nCols=2;
hf=figure;
fpos=get(hf,'Position');
set(hf,'Position',[fpos(1) fpos(2) fpos(3) fpos(4)]);
hs(1)=subplot(nRows,nCols,1); hold on
for e=1:3
    plot(time,grandMeanPos(:,e),'LineWidth',1.5);
    %hsh(e,1)=shadedErrorBar(time,grandMuPos(:,e),grandSemPos(:,e),{'markerfacecolor',colors(e,:)},1);
end
axis tight
set(get(hs(1),'Xlabel'),'String','Time (mm:ss)');
set(hs(1),'Xtick',[0 600 1200 1800]);
set(hs(1),'XtickLabel',{'0:00','10:00','20:00','30:00'});
yl=ylim;
set(hs(1),'YTick',-0.2:0.2:0.2);
plot([600 600],[yl(1) yl(2)],'-k');
plot([1200 1200],[yl(1) yl(2)],'-k');
harea=area([600 1200],[yl(2) yl(2)],'BaseValue',yl(1),'FaceColor',[0.7 0.7 0.7],'FaceAlpha',0.15,'LineStyle','none');
%htit(5)=title('\beta>0','FontWeight','normal');
set(get(hs(1),'ylabel'),'String','BOLD (a.u.)');
htext3=text(30,0.25,'$\beta>0$');set(htext3,'Interpreter','Latex');set(htext3,'FontSize',14);
htext3=text(690,-0.2,'Light on'); set(htext3,'FontSize',14);
hlg=legend('Echo 1','Echo 2','Echo 3');
set(hlg,'Box','off');
set(hlg,'orientation','horizontal');
lgPos=get(hlg,'Position');
set(hlg,'Position',[lgPos(1)+0.025 lgPos(2)+0.05 lgPos(3) lgPos(4)]);

hs(2)=subplot(nRows,nCols,2); hold on
for e=1:3
    plot(time,grandMeanNeg(:,e),'LineWidth',1.5);
    %hsh(e,2)=shadedErrorBar(time,grandMuNeg(:,e),grandSemNeg(:,e),{'markerfacecolor',colors(e,:)},1);
end
axis tight
set(get(hs(2),'Xlabel'),'String','Time (mm:ss)');
set(hs(2),'Xtick',[0 600 1200 1800]);
set(hs(2),'XtickLabel',{'0:00','10:00','20:00','30:00'});
set(hs(2),'ylim',yl);
set(hs(2),'YTick',-0.2:0.2:0.2);
set(hs(2),'YTickLabel',{'','',''});
plot([600 600],[yl(1) yl(2)],'-k');
plot([1200 1200],[yl(1) yl(2)],'-k');
harea=area([600 1200],[yl(2) yl(2)],'BaseValue',yl(1),'FaceColor',[0.7 0.7 0.7],'FaceAlpha',0.15,'LineStyle','none');
%htit(5)=title('\beta<0','FontWeight','normal');
htext4=text(30,0.25,'$\beta<0$');set(htext4,'Interpreter','Latex');set(htext4,'FontSize',14);
htext4=text(690,0.22,'Light on'); set(htext4,'FontSize',14);




%% 

hs(3)=subplot(nRows,nCols,3); hold on
for e=1:3
    hbar(e,1)=bar(e,mubpos(e));
    herr=errorbar([1 2 3],mubpos,sembpos,'LineStyle','none','Color','k');
    for e=1:3
        yl=ylim;
        if pvalpos(e)<0.05
            text(e-0.01, yl(1) + (yl(2)-yl(2))*0.9, '*');
        end
    end
end
set(gca,'XAxisLocation','top');

set(gca,'Xtick',[1 2 3]);
set(gca,'Ytick',[-0.001 0]);
set(gca,'YTickLabel',{'-0.001','0'});
hyl=ylabel('Baseline slope (units/s)');
%xlabel('Echo');
set(gca,'XtickLabel',{'','',''});
htext3=text(-0.1,0.0001,'$\beta>0$');set(htext3,'Interpreter','Latex');set(htext3,'FontSize',14);
%xlim([0.5 3.5]);
ylpos=get(hyl,'Position');
set(hyl,'Position',[-1.5 -0.00115 ylpos(3)]);
box off

hs(4)=subplot(nRows,nCols,5); hold on
for e=1:3
    hbar(e,1)=bar(e,mubneg(e)); 
    herr=errorbar([1 2 3],mubneg,sembneg,'LineStyle','none','Color','k');
    for e=1:3
        yl=ylim;
        if pvalneg(e)<0.05
            text(e-0.01, yl(1) + (yl(2)-yl(1))*1.05, '*');
        end
    end
end

set(gca,'Xtick',[1 2 3]);
set(gca,'Ytick',[0 0.001]);
set(gca,'YTickLabel',{'0','0.001'});
xlabel('Echo');
htext3=text(-0.1,0.0011,'$\beta<0$');set(htext3,'Interpreter','Latex');set(htext3,'FontSize',14);
%xlim([0.5 3.5]);
%%
allposx=[]; allposy=[];
allnegx=[]; allnegy=[];
allpose=[]; allnege=[];
% scatter of slope vs effect size?
for e=1:3
    tx=bpos(e,~isnan(bpos(e,:)));
    ty=mupscpos(e,~isnan(mupscpos(e,:)));
    allposx=cat(1,allposx,tx(:));
    allposy=cat(1,allposy,ty(:));
    allpose=cat(1,allpose,e*ones(numel(tx(:)),1));
    
    tx=bneg(e,~isnan(bneg(e,:)));
    ty=mupscneg(e,~isnan(mupscneg(e,:)));
    allnegx=cat(1,allnegx,tx(:));
    allnegy=cat(1,allnegy,ty(:));
    allnege=cat(1,allnege,e*ones(numel(tx(:)),1));
end

[r,ppos]=corrcoef(allposx,allposy);
[r,pneg]=corrcoef(allnegx,allnegy);
[rall,pall]=corrcoef([allposx; allnegx],[allposy; allnegy]);

%%
hs(5)=subplot(nRows,nCols,nRows*nCols); hold on
%hsc=scatter([allposx; allnegx],[allposy; allnegy]);
for e=1:3
    % positive scatters
    hsc(e,1)=scatter(allposx(allpose==e),allposy(allpose==e),'MarkerFaceColor',colors(e,:),'Marker','^','MarkerEdgeColor','none');
    % negative scatters
    hsc(e,2)=scatter(allnegx(allnege==e),allnegy(allnege==e),'MarkerFaceColor',colors(e,:),'Marker','v','MarkerEdgeColor','none');
end
ylim([-0.2 0.2]);
xlim([-0.002 0.002]);
xl=xlim;
yl=ylim;
set(gca,'Xtick',[-0.002 0 0.002]);
set(gca,'Ytick',[-0.2 0 0.2]);
set(gca,'XTickLabel',{'-0.002','0','0.002'});
Y=[allposy; allnegy];
X=[allposx; allnegx];
Bfit=regress(Y,[X,ones(size(X,1),1)]);
hfit=plot([xl(1) xl(2)],[xl(1) xl(2)]*Bfit(1)+Bfit(2),'--k');
%htxt=text(0.001,0.1,sprintf(['r=%0.2f , p=%0.2g'],rall(1,2),pall(1,2)));
htxt=text(0.0005,0.1,sprintf(['r=%0.2f, p<0.001'],rall(1,2)));
ylabel('Beta weight');
hxl=xlabel('Baseline slope (units/second)');

%%
spos=get(hs(1),'Position');
set(hs(1),'Position',[spos(1)-0.025 spos(2) spos(3)*1.2 spos(4)]);
spos=get(hs(2),'Position');
set(hs(2),'Position',[spos(1) spos(2) spos(3)*1.2 spos(4)]);

spos=get(hs(3),'Position');
set(hs(3),'Position',[spos(1)+0.05 spos(2)+0.06 spos(3)*0.7 spos(4)*0.6]);
spos=get(hs(4),'Position');
set(hs(4),'Position',[spos(1)+0.05 spos(2)+0.175 spos(3)*0.7 spos(4)*0.6]);
% 
spos=get(hs(5),'Position');
set(hs(5),'Position',[spos(1) spos(2)+0.175 spos(3)*1.15 spos(4)*1.45]);

%%
sublabel([hs(1) hs(2) hs(3) hs(5)],-10,-33,'FontWeight','Bold','FontSize',16);

if PRINT_FIGURE
    %%
    %
    %  PREFORMATTED
    %  TEXT
    %
    figStr=['../figures/revfig1b-' date '.png'];
    print('-dpng',figStr); % ../figures/revfig1
    crop(figStr,0);
end