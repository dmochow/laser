%% 05.21.19
% show the time courses of all subjects

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


load(savedDataFilename,'r','r2','allts','issig','allBetas','allOnsets');

subj2show=find(sum(r)>0);
nSubjShow=numel(subj2show);
nRows=nSubjShow/2; 
allsigposts=cell(3,1);
allsignegts=cell(3,1);
%%
figure
for s=1:nSubjShow
    onsetTR=floor(allOnsets(subj2show(s),1)/TR);
    tshift=215-onsetTR;
    hs(s)=subplot(nRows,2,s); hold on
    for e=1:3
        allts{e,s}=circshift(allts{e,s},[tshift 0]);
        allsigposts{e}=cat(2,allsigposts{e},allts{e,s}(:,issig{e,s} & allBetas{e,s}(1,:)>0));
        allsignegts{e}=cat(2,allsignegts{e},allts{e,s}(:,issig{e,s} & allBetas{e,s}(1,:)<0));
        
        sigposts{e,s}=allts{e,s}(:,issig{e,s} & allBetas{e,s}(1,:)>0);
        signegts{e,s}=allts{e,s}(:,issig{e,s} & allBetas{e,s}(1,:)<0);
        muposts(:,e,s)=mean(sigposts{e,s},2);
        munegts(:,e,s)=mean(signegts{e,s},2);
        
%         allts{e,subj2show(s)}=circshift(allts{e,subj2show(s)},[tshift 0]);
%         thisTs=allts{e,subj2show(s)};
%         vinds=issig{e,subj2show(s)};
%         muThisTs=mean(thisTs(:,vinds==1),2);
        plot(time,muposts(:,e,s));
        
        
        
    end
    %ylim([-1 1]);
    axis tight
    yl=ylim;

    plot([600 600],[yl(1) yl(2)],'-k');
    plot([1200 1200],[yl(1) yl(2)],'-k');
    harea=area([600 1200],[yl(2) yl(2)],'BaseValue',yl(1),'FaceColor',[0.7 0.7 0.7],'FaceAlpha',0.15,'LineStyle','none');

    htxt=text(50,yl(1)+(yl(2)-yl(1))*0.1,['S' num2str(subj2show(s))]);
    
    if s==1, 
        hlg=legend('Echo 1','Echo 2','Echo 3'); 
        htxt=text(725,yl(1)+0.125,['Light on']);
    end
    
end


%%
% esthetics

% 
set(hs(1:end),'XLim',[0 1800]);
set(hs(1:end),'XTick',[0 600 1200 1800]);
set(hs(1:end),'XTickLabel',{'','','',''});
set(hs(end-1),'XTickLabel',{'0','10','20','30'});
set(hs(end),'XTickLabel',{'0','10','20','30'});
% 
set(hs(1:end), 'YLimSpec', 'Tight');

set(hs(1:end),'Ytick',[-0.25 0 0.25]);
set(hs(2:2:end),'YtickLabel',{'','',''});
%set(hs(3:2:end),'YTickLabel',{'','',''});
% 
set(get(hs(1),'YLabel'),'String','BOLD (a.u.)');
set(get(hs(end-1),'XLabel'),'String','Time (min)');
set(get(hs(end),'XLabel'),'String','Time (min)');

% 
ax=1.2; ay=1.2;
for s=1:nSubjShow
    pos=get(hs(s),'Position');
    set(hs(s),'Position',[pos(1) pos(2) pos(3)*ax pos(4)*ay]);
end
% 
lgPos=get(hlg,'Position');
delx=-0.12; dely=0.07;
set(hlg,'Position',[lgPos(1)+delx lgPos(2)+dely lgPos(3) lgPos(4)]);
set(hlg,'box','off','Orientation','horizontal');
% 
% 
if PRINT_FIGURE
    %%
    % 
    %  PREFORMATTED
    %  TEXT
    % 
    figStr=['../figures/allTimeCourses.png'];
    print('-dpng',figStr);
    crop(figStr,0);
end

