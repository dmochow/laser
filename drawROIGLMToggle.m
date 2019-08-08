%% 06.01.19
% draw results of Study II (toggling)
clear all; close all; clc
PRINT_FIGURE=1;
savedDataFilename='../data/precomputed/revFig2data.mat';
subjStrs={'S23','S24','S25','S26','S27','S28','S29','S30','S31','S32'};
nSubjects=numel(subjStrs);
TR=2.8;
alpha=0.05; % significance level
gmIndex=2; % this will never change

load(savedDataFilename,'issig','allBetas','allts','allRegressors','r','r2');

%%
% construct grand mean time course for significant voxels
% split based on positive or negative beta
allsigposts=cell(3,1);
allsignegts=cell(3,1);
for s=1:nSubjects
    for e=1:3        
        % figure out positive or negative betas
        % note: using B(1) here (acute effect)
        % B(2) is constant term
        tbetas=allBetas{e,s}(1,:);
        tbetas(~issig{e,s})=0; % remove non-significant ones
        propSigPos(e,s)= mean(tbetas>0);
        propSigNeg(e,s)= mean(tbetas<0);
        sumSigPos(e,s)= numel(find(tbetas>0));
        sumSigNeg(e,s)= numel(find(tbetas<0));
        
        sigbetaspos3=allBetas{e,s}(:,find(issig{e,s} & allBetas{e,s}(1,:)>0   ));
        pscpos{e,s}=100*abs(sigbetaspos3(1,:)./sigbetaspos3(2,:));
        pscposb{e,s}=sigbetaspos3(1,:);
        
        sigbetasneg3=allBetas{e,s}(:,find(issig{e,s} & allBetas{e,s}(1,:)<0   ));
        pscneg{e,s}=100*abs(sigbetasneg3(1,:)./sigbetasneg3(2,:));
        pscnegb{e,s}=sigbetasneg3(1,:);
        
    end
end

%% # SIGNIFICANT
figure
hs(1)=subplot(3,1,1); hold on
hbar=bar(1:10,[100*propSigPos' -100*propSigNeg']);
xl=xlim; 
ylim([-40 40]);
yl=ylim;
harea1=area([xl(1) xl(2)],[yl(2) yl(2)],'BaseValue',0,'FaceColor',[0.6 0.6 0.6],'FaceAlpha',0.15,'LineStyle','none');
harea2=area([xl(1) xl(2)],[yl(1) yl(1)],'BaseValue',0,'FaceColor',[0.3 0.3 0.3],'FaceAlpha',0.15,'LineStyle','none');
htext1=text(1,30,'$\beta>0$');set(htext1,'Interpreter','Latex');set(htext1,'FontSize',14);
htext2=text(1,-30,'$\beta<0$');set(htext2,'Interpreter','Latex');set(htext2,'FontSize',14);
set(get(hs(1),'ylabel'),'String','% sig voxels');
set(hs(1),'Xtick',1:10);
set(hs(1),'XTickLabel',{'','','','','', ...
    '','','','',''});
%set(hs(1),'XTickLabelRotation',45);
set(hs(1),'Ytick',[-40 0 40]);
set(hs(1),'YTickLabel',{'40','0','40'});
hlg0=legend('Echo 1','Echo 2','Echo 3');
set(hlg0,'Orientation','horizontal');
set(hlg0,'box','off');
set(hbar(4),'FaceColor',get(hbar(1),'FaceColor'));
set(hbar(5),'FaceColor',get(hbar(2),'FaceColor'));
set(hbar(6),'FaceColor',get(hbar(3),'FaceColor'));
ylim(yl);

%% EFFECT SIZE
muBetasPos=cellfun(@mean,pscposb)*1;
sigmaBetasPos=cellfun(@std,pscposb)*1;
muBetasNeg=cellfun(@mean,pscnegb)*1;
sigmaBetasNeg=cellfun(@std,pscnegb)*1;
hs(2)=subplot(3,1,2); hold on
hbar=bar(1:10,[muBetasPos' muBetasNeg']);

ngroups=10;
nbars=6;
groupwidth = min(0.8, nbars/(nbars + 1.5));
y=[muBetasPos' muBetasNeg'];
err=[sigmaBetasPos' sigmaBetasNeg'];
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    herr(:,i)=errorbar(x,y(:,i), err(:,i), '.k','CapSize',0);
end
%hold off

xl=xlim; yl=ylim;
harea1=area([xl(1) xl(2)],[yl(2) yl(2)],'BaseValue',0,'FaceColor',[0.6 0.6 0.6],'FaceAlpha',0.15,'LineStyle','none');
harea2=area([xl(1) xl(2)],[yl(1) yl(1)],'BaseValue',0,'FaceColor',[0.3 0.3 0.3],'FaceAlpha',0.15,'LineStyle','none');

xl=get(hs(1),'xlim'); set(hs(2),'XLim',xl);
set(get(hs(2),'Ylabel'),'String','Beta weights');
set(hs(2),'Xtick',1:20);
set(hs(2),'XTickLabel',{'S1','S2','S3','S4','S5', ...
    'S6','S7','S8','S9','S10'});
set(hs(2),'XTickLabelRotation',45);
set(hs(2),'Ytick',[-0.2 0 0.2]);

set(hbar(4),'FaceColor',get(hbar(1),'FaceColor'));
set(hbar(5),'FaceColor',get(hbar(2),'FaceColor'));
set(hbar(6),'FaceColor',get(hbar(3),'FaceColor'));

%%
%look at the significant voxels
duration=round(180/TR);
cycleTime=(0:duration-1)*TR;
for e=1:3
    for s=1:nSubjects
                regressor=allRegressors{s};
                ontimes=find(diff(regressor)==1);
                
        sigvoxNeg=find(issig{e,s} & allBetas{e,s}(1,:)<0 );
        tcsNeg=allts{e,s}(:,sigvoxNeg);
  
        nMissing=ontimes(end)+duration-size(tcsNeg,1);
        tcsNeg=cat(1,tcsNeg,zeros(nMissing,size(tcsNeg,2)));
        epochsNeg{e,s} = mean ( simpleEpoch(tcsNeg',ontimes,duration) , 3); % mean across cycles
        muEpochsNeg(:,e,s) = mean ( epochsNeg{e,s} ); % mean across voxels
        
        sigvoxPos=find(issig{e,s} & allBetas{e,s}(1,:)>0 );
        tcsPos=allts{e,s}(:,sigvoxPos);
        nMissing=ontimes(end)+duration-size(tcsPos,1);
        tcsPos=cat(1,tcsPos,zeros(nMissing,size(tcsPos,2)));
        epochsPos{e,s} = mean ( simpleEpoch(tcsPos',ontimes,duration) , 3); % mean across cycles
        muEpochsPos(:,e,s) = mean ( epochsPos{e,s} ); % mean across voxels
    end
end

%%
[semCycleNeg,muCycleNeg]=nansem(muEpochsNeg,3);
[semCyclePos,muCyclePos]=nansem(muEpochsPos,3);

hs(3)=subplot(3,2,5); hold on
for e=1:3
    hsh(e)=shadedErrorBar(cycleTime,muCyclePos(:,e),semCyclePos(:,e),{'markerfacecolor',get(hbar(e),'FaceColor')},1);
end
yl=ylim;
plot([0 0],[yl(1) yl(2)],'-k');
plot([60 60],[yl(1) yl(2)],'-k');
harea=area([0 60],[yl(2) yl(2)],'BaseValue',yl(1),'FaceColor',[0.7 0.7 0.7],'FaceAlpha',0.15,'LineStyle','none');
htit1=title('Cycle Average (\beta>0)','FontWeight','normal');
%hlg3=legend(harea,'Light on');
%set(hlg3,'box','off');
set(hs(3),'xtick',[0 60 120 180]);
set(hs(3),'XTickLabel',{'0:00','1:00','2:00','3:00'});
set(hs(3),'ytick',[-0.1 0 0.1]);
set(get(hs(3),'ylabel'),'String','BOLD (a.u.)');
htext3=text(3,-0.09,'Light on'); set(htext3,'FontSize',14);
set(get(hs(3),'Xlabel'),'String','Time (s)');

hs(4)=subplot(3,2,6); hold on
for e=1:3
    hsh(e)=shadedErrorBar(cycleTime,muCycleNeg(:,e),semCycleNeg(:,e),{'markerfacecolor',get(hbar(e),'FaceColor')},1);
end
ylim(yl);
plot([0 0],[yl(1) yl(2)],'-k');
plot([60 60],[yl(1) yl(2)],'-k');
harea=area([0 60],[yl(2) yl(2)],'BaseValue',yl(1),'FaceColor',[0.7 0.7 0.7],'FaceAlpha',0.15,'LineStyle','none');
htit1=title('Cycle Average (\beta<0)','FontWeight','normal');
set(hs(4),'ytick',[-0.1 0 0.1]);
htext4=text(3,0.1,'Light on'); set(htext4,'FontSize',14);
set(get(hs(4),'Xlabel'),'String','Time (s)');

% hins = axes('position', [0.15 0.91 0.25 0.05]);
% plot(hins,[zeros(1,666) repmat( [ones(1,333) zeros(1,666)], 1 ,  10   )],'k','LineWidth',1.5); box off
% htxt=text(360,1.37,'Design');
% set(hins,'XTick',[]); set(hins,'YTick',[]);

subpos=get(hs(1),'Position');
set(hs(1),'Position',[subpos(1) subpos(2) subpos(3) subpos(4)]);

subpos=get(hs(2),'Position');
set(hs(2),'Position',[subpos(1) subpos(2)+0.05 subpos(3) subpos(4)]);

subpos=get(hs(3),'Position');
set(hs(3),'Position',[subpos(1) subpos(2) subpos(3) subpos(4)]);

subpos=get(hs(4),'Position');
set(hs(4),'Position',[subpos(1) subpos(2) subpos(3) subpos(4)]);

%% move subplots 1 and 2 to make room for design matrix
for s=1:2
    spos=get(hs(s),'Position');
    set(hs(s),'Position',[0.33 spos(2) spos(3)*0.8 spos(4)]);
end

%freezeColors
%% DESIGN MATRIX
hs(5)=subplot(251);  hold on
designMatrix2(:,1)=[zeros(1,333) repmat([ones(1,333) zeros(1,666)],1,10)];
designMatrix2(:,2)=[zeros(1,333) repmat([ones(1,333) zeros(1,666)],1,10)];
himg=imagesc([1 2],[0 1800],designMatrix2); 
cm=colormap; 
cm(1,:)=[0 0 0];
cm(end,:)=[1 0 0];
colormap(cm([1 end],:));
axis ij;
htits=title('Design Matrix','FontWeight','normal'); 
set(hs(5),'Ytick',[0 600 1200 1800]); 
set(hs(5),'YTickLabel',{'0','10','20','30'}); 
set(hs(5),'XTick',[1.25 1.75]);set(hs(3),'XTickLabel',{'1','2'}); 
set(hs(5),'XTick',[]);
xl=xlim; yl=ylim;
xlim([1 2]); ylim([0 1800]);
ylabel('Time (min)');
subpos=get(hs(5),'Position');
set(hs(5),'Position',[subpos(1)-0.025 subpos(2)-0.08 subpos(3)*1.25 subpos(4)*1.25]);
hcb=colorbar; set(hcb,'Ticks',[0.25 0.75]); set(hcb,'TickLabels',{'0','1'});
titpos=get(htits,'Position');
set(htits,'Position',[1.5000 2000 0]);
%axis square
%%

if PRINT_FIGURE
    %%
    %
    %  PREFORMATTED
    %  TEXT
    %
    sublabel([hs(5) hs(1) hs(2)  hs(3) hs(4)],0,-40,'FontWeight','Bold','FontSize',16);
    figStr=['../figures/revfig2-' date '.png'];
    print('-dpng',figStr); % ../figures/revfig1
    crop(figStr,0);
end


%% stats for reporting in paper
x=abs([muBetasPos muBetasNeg]);
mu=nanmean(x,2)
sigma=nanstd(x,[],2)


