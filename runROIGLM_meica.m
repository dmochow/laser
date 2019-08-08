%% 05.04.19
% MEICA preprocessing
% analysis script for Study I
% GLM analysis of BOLD near the illumination site
clear all; close all; clc
PRINT_FIGURE=0;

subjStrs={'S02','S03','S04','S05','S06','S07','S09','S10','S11','S12','S13'...
    ,'S14','S15','S16','S17','S18','S19','S20','S21','S22'};
nSubjects=numel(subjStrs);

subjectNames={'Greg','Emircan','Pat','Flores','Christian',...
    'Jamal','Kevin','Destiny','Jacek','Marcelo',...
    'James','Enisena','Felicia','Ridmila','Karem',...
    'Blerita','Duc','Anu','Javier','Henintsoa'};

TR=2.8;
nTRs=643; % hard-coded but needs to be changed if this script is to be adapted for the toggling study
alpha=0.05; % significance level
gmIndex=2; % this will never change
roiMaskFilename='follow_ROI_LASER+tlrc.BRIK';
altRoiMaskFilename='mu_follow_ROI_LASER+tlrc.BRIK'; % no markers
segFilename='resampled_Classes+tlrc.BRIK';
%boldPrefix='all_runs.S02+tlrc.BRIK'; % use this to test different preprocs
biopacFilename='biopac.mat';

for s=1:nSubjects
    s
    pathToData=['../data/output/' subjStrs{s} '/' subjStrs{s} '.results'];
    pathToBiopac=['../data/scanner_data/' subjStrs{s}];
    
    boldFilename=['all_runs.' subjStrs{s} '+tlrc.BRIK']; % use this to test different preprocs
    
    % illumination ROI
    [~, roiMask, info, ~] = BrikLoad (fullfile(pathToData,roiMaskFilename));
    roiMask=logical(roiMask);
    
    if isempty(roiMask) % not all subjects had markers
        [~, roiMask, info, ~] = BrikLoad (fullfile(pathToData,altRoiMaskFilename));
        roiMask=logical(roiMask);
    end
    
    % grey matter mask
    resampledSegFilename=fullfile(pathToData,'resampled_Classes+tlrc.BRIK');
    [~, seg, ~, ~] = BrikLoad (resampledSegFilename);
    greyMask=logical(seg==gmIndex);
    
    finalMask=roiMask&greyMask;
    
    % TODO: load biopac here
    % define onset and offset time
    load(fullfile(pathToBiopac,biopacFilename),'data');
    [~,maxindsample]=max(diff(data(:,1)));
    onsetTimeSec=maxindsample/1000;
    onsetTimeTR=round(onsetTimeSec/TR);
    [~,minindsample]=min(diff(data(:,1)));
    offsetTimeSec=minindsample/1000;
    offsetTimeTR=round(offsetTimeSec/TR);
    allOnsets(s,:)=[onsetTimeTR,offsetTimeTR]*TR;
    
    t2=[onsetTimeTR:offsetTimeTR-1];
    t3=[offsetTimeTR:nTRs];
    X=zeros(nTRs,2);
    X(t2,1)=1;
    X(t3,2)=1;
    Xnull=ones(nTRs,1);
    
    
    [~, bold, boldInfo, ~] = BrikLoad (fullfile(pathToData,boldFilename));
    
    %%
    ts=vol2ts(bold,finalMask);
    allts{s}=ts;
    nVoxels=size(ts,2);
    for v=1:nVoxels
        statsf=myGLM(ts(:,v),X,1); % full model
        statsr=myGLM(ts(:,v),Xnull,0); % reduced model
        allBetas{s}(v)=statsf.B(1);
        ssef=statsf.sse;
        doff=statsf.dof;
        sser=statsr.sse;
        dofr=statsr.dof;
        Fstat=((sser-ssef)/(dofr-doff)) / (ssef/doff);
        pval{s}(v)=fcdf(Fstat,dofr-doff,doff,'upper');
    end
    [p_fdr,p_masked]=fdr(pval{s},alpha);
    issig{s}=p_masked;
    
end



%%
r=cellfun(@mean,issig);
r2=cellfun(@sum,issig);

%%
figure
hs(1)=subplot(311);
bar(1:20,100*r');
set(get(hs(1),'ylabel'),'String','% significant voxels');
set(get(hs(1),'xlabel'),'String','Subject');
hlg0=legend('Echo 1','Echo 2','Echo 3');

%
%figure
% find best subjects
[~,sortind]=sort(r,'descend');
nSubjShow=8;
TR=2.8;
time=(0:642)*TR;
for s=1:nSubjShow
    hs(s+1)=subplot(6,2,4+s); hold on
    thisTs=allts{sortind(s)};
    vinds=issig{sortind(s)};
    muThisTs=mean(thisTs(:,vinds==1),2);
    plot(time,muThisTs);
    %ylim([-1 1]);
    axis tight
    yl=ylim;
    plot([allOnsets(s,1) allOnsets(s,1)],[yl(1) yl(2)],'-k');
    plot([allOnsets(s,2) allOnsets(s,2)],[yl(1) yl(2)],'-k');
    harea=area([allOnsets(s,1) allOnsets(s,2)],[yl(2) yl(2)],'BaseValue',yl(1),'FaceColor',[0.7 0.7 0.7],'FaceAlpha',0.15,'LineStyle','none');
    %plot([600 600],[yl(1) yl(2)],'-k');
    %plot([1200 1200],[yl(1) yl(2)],'-k');
    %harea=area([600 1200],[yl(2) yl(2)],'BaseValue',yl(1),'FaceColor',[0.7 0.7 0.7],'FaceAlpha',0.15,'LineStyle','none');
    
    %htxt=text(50,yl(2)-0.25,['S' num2str(sortind(s))]);
    htxt=text(50,yl(1)+0.15,['S' num2str(sortind(s))]);
    
%     if s==1,
%         hlg=legend('Echo 1','Echo 2','Echo 3');
%         htxt=text(725,yl(1)+0.25,['Light on']);
%     end
    
end

hins = axes('position', [0.14 0.85 0.2 0.05]);
plot(hins,1:999,[zeros(1,333) ones(1,333) zeros(1,333)],'k','LineWidth',2); box off
htxt=text(360,1.37,'Design');
set(hins,'XTick',[]); set(hins,'YTick',[]);

% esthetics
set(hs(1),'Xtick',1:20);
set(hs(1),'Box','off');

set(hs(2:end),'XLim',[0 1800]);
set(hs(2:end),'XTick',[0 600 1200 1800]);
set(hs(2:nSubjShow-1),'XTickLabel',{'','','',''});
set(hs(nSubjShow),'XTickLabel',{'0:00','10:00','20:00','30:00'});
set(hs(nSubjShow+1),'XTickLabel',{'0:00','10:00','20:00','30:00'});

set(hs(2:end), 'YLimSpec', 'Tight');

%set(hs(2:end),'Ytick',[-1 0 1]);
%set(hs(3:2:end),'YTickLabel',{'','',''});

set(get(hs(2),'YLabel'),'String','BOLD (a.u.)');
set(get(hs(nSubjShow),'XLabel'),'String','Time (mm:ss)');

set(hlg0,'box','off','Orientation','horizontal');

lgPos=get(hlg,'Position');
delx=-0.07; dely=-0.07;
set(hlg,'Position',[lgPos(1)+delx lgPos(2)+dely lgPos(3) lgPos(4)]);
set(hlg,'box','off','Orientation','horizontal');

sublabel([hs(1) hs(2)],-10,-60,'FontWeight','Bold','FontSize',16);

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
%
%   for x = 1:10
%       disp(x)
%   end
%

%% %%
% s=7;
% b1=allBetas{s}(issig{s});
% b2=allBetas{s}(issig{s});
% b3=allBetas{s}(issig{s});
% figure;
% hold on
% plot(ones(1,numel(b1)),b1,'o');
% plot(2*ones(1,numel(b2)),b2,'o');
% plot(3*ones(1,numel(b3)),b3,'o');