clear all; close all; clc
subjStrs={'S23','S24','S25','S26','S27','S28','S29','S30','S31','S32'};
nSubjects=numel(subjStrs);
alpha=0.05;
nTRs=685; % nTRs in the longest recording
TR=2.8;
anatRoiMaskFilename='roi_r21_z39+tlrc.BRIK';
anatBrainMaskFilename='brain_mask+tlrc.BRIK';
precomputedFilename='../data/precomputed/ROIGLMToggleStats.mat';


%%
% load data
load(precomputedFilename,'issig','allBetas','allts','allRegressors','r','r2');
r=r(:,1:nSubjects); % hack as we collect data

%%
% figure
hf=figure;
nRows=4;
% %significant vs subject
for s=1:nSubjects
    
    pathToData=['../data/myoutput/' subjStrs{s}];
    [~, anat, ~, ~] = BrikLoad (fullfile(pathToData,'anat+tlrc'));
    anat=anat(end:-1:1,:,:);
    
    [~, anatRoiMask, ~, ~] = BrikLoad (fullfile(pathToData,anatRoiMaskFilename));
    anatRoiMask=anatRoiMask(end:-1:1,:,:);
    anatRoiMask=logical(anatRoiMask);
    
    %zcoor=75;
    stats = regionprops(anatRoiMask);
    zcoor = round(stats.Centroid(3));
    
    [~, anatBrainMask, ~, ~] = BrikLoad (fullfile(pathToData,anatBrainMaskFilename));
    anatBrainMask=anatBrainMask(end:-1:1,:,:);
    
    for e=1:3
        
        [~, issigMask, ~, ~] = BrikLoad (fullfile(pathToData,['issigAnatEcho' num2str(e) '+tlrc']));
        issigMask=issigMask(end:-1:1,:,:);
        
        %hs(1,s)=subplot(nRows,nSubjects,(e-1)*nSubjects+s);
        hs(e,s)=subplot(nRows,nSubjects,e*nSubjects+s);
        h = colorBrain(anat(:,:,zcoor),anatRoiMask(:,:,zcoor),issigMask(:,:,zcoor)==1,anatBrainMask(:,:,zcoor),1);
        set(h,'FaceAlpha',0.9);
    end
    
    
    
end

%hs(nRows,1)=subplot(nRows,1,nRows);
hs(nRows,1)=subplot(nRows,1,1);
hbar=bar(1:nSubjects,100*r');
box off
set(get(hs(nRows,1),'ylabel'),'String','% significant voxels');
set(get(hs(nRows,1),'xlabel'),'String','Subject');
hlg0=legend('Echo 1','Echo 2','Echo 3'); 
set(hlg0,'box','off','Orientation','horizontal');
lgPos=get(hlg0,'Position');
delx=0; dely=0.04;
set(hlg0,'Position',[lgPos(1)+delx lgPos(2)+dely lgPos(3) lgPos(4)]);


hins = axes('position', [0.15 0.89 0.25 0.05]);
plot(hins,[zeros(1,666) repmat( [ones(1,333) zeros(1,666)], 1 ,  10   )],'k','LineWidth',1.5); box off
htxt=text(360,1.37,'Design');
set(hins,'XTick',[]); set(hins,'YTick',[]);

%%
dely=-0.05;
for s=1:nSubjects
    pos=get(hs(2,s),'Position');
    set(hs(2,s),'Position',[pos(1) pos(2)+dely pos(3) pos(4)]);
end

%%
dely=-0.1;
for s=1:nSubjects
    pos=get(hs(1,s),'Position');
    set(hs(1,s),'Position',[pos(1) pos(2)+dely pos(3) pos(4)]);
end

htxt1=text(-3200,-7,'Echo 1');
htxt2=text(-3200,-10.5,'Echo 2');
htxt3=text(-3200,-14,'Echo 3');

sublabel([hs(4,1) hs(1,1)],-10,-60,'FontWeight','Bold','FontSize',16);

print -dpng ../figures/revfig2
crop('../figures/revfig2.png',0);

%% time course figure
TR=2.8;
duration=round(180/TR)-2;
figure;
nRows=nSubjects; nCols=2;
for s=1:nSubjects
    regressor=allRegressors{s};
    ontimes=find(diff(regressor)==1);
    
    hs(s,1)=subplot(nRows,nCols,(s-1)*2+1); hold on
    for e=1:3
        sigvox=find(issig{e,s});
        tcs=allts{e,s}(:,sigvox);
        time=(0:size(tcs,1)-1)*TR;
        plot(time,mean(tcs,2));
        yl=ylim;
        for o=1:numel(ontimes)
            tstart=ontimes(o)*TR; tend=(ontimes(o)+21)*TR;
            harea=area([tstart tend],[yl(2) yl(2)],'BaseValue',yl(1),'FaceColor',[0.7 0.7 0.7],'FaceAlpha',0.15,'LineStyle','none');
        end
        %plot(regressor*0.5,'k');
    end
    
    hs(s,2)=subplot(nRows,nCols,s*2); hold on
    time2=(0:61)*TR;
    for e=1:3
        sigvox=find(issig{e,s});
        tcs=allts{e,s}(:,sigvox);
        ontimes=find(diff(regressor)==1);
        epochs = simpleEpoch(tcs',ontimes,duration);
        muEpoch=mean(epochs,3);
        muMuEpoch=mean(muEpoch,1);
        plot(time2,muMuEpoch);
        yl=ylim;
        harea=area([0 60],[yl(2) yl(2)],'BaseValue',yl(1),'FaceColor',[0.7 0.7 0.7],'FaceAlpha',0.15,'LineStyle','none');
    end
    
end

set(hs(:,1),'XLim',[0 1800]);
set(hs(:,1),'XTick',[0 600 1200 1800]);
set(hs(1:end-1,1),'XTickLabel',{'','','',''});
set(hs(end,1),'XTickLabel',{'0:00','10:00','20:00','30:00'});
%set(hs(nSubjShow+1),'XTickLabel',{'0:00','10:00','20:00','30:00'});


set(hs(:,2),'XLim',[0 180]);
set(hs(:,2),'XTick',[0 60 120 180]);
set(hs(1:end-1,2),'XTickLabel',{'','','',''});
set(hs(end,2),'XTickLabel',{'0:00','1:00','2:00','3:00'});



for s=1:nSubjects
    set(get(hs(s,1),'Title'),'String',['Subject ' num2str(s) ]);
    set(get(hs(s,1),'Title'),'FontWeight','normal');
end

set(get(hs(1,2),'Title'),'String','Cycle Average');
set(get(hs(1,2),'Title'),'FontWeight','normal');

print -dpng ../figures/tmp

%%% JUNK
%%
% % look at the significant voxels
% duration=round(180/TR)-2;
% figure; hold on
% for e=1:3
%     for s=1:nSubjects
%         sigvox=find(issig{e,s});
%         tcs=allts{e,s}(:,sigvox);
%         %plot(mean(tcs,2)); plot(regressor);
%         ontimes=find(diff(regressor)==1);
%         epochs = simpleEpoch(tcs',ontimes,duration);
%         muEpoch=mean(epochs,3);
%         muMuEpoch=mean(muEpoch,1);
%         hs(e,s)=subplot(3,nSubjects,(e-1)*nSubjects+s);
%         plot(muMuEpoch);
%     end
% end

% %%
% allsigts=cell(3,1);
% for e=1:3
%     for s=1:nSubjects
%         tmp=allts{e,s}(:,issig{e,s});
%         nn=size(tmp,1);
%         if nn < nTRs
%             tmp=cat(1,tmp,zeros(nTRs-nn,size(tmp,2)));
%         end
%         allsigts{e}=cat(2,allsigts{e},tmp);
%     end
% end
%
% %%
% figure;
% for e=1:3
%     subplot(3,1,e);
%     plot(mean(allsigts{e},2));
% end

% % two time courses
% s=4;
% hs(2)=subplot(5,2,3); hold on
% for e=1:3
%     sigvox=find(issig{e,s});
%     tcs=allts{e,s}(:,sigvox);
%     plot(mean(tcs,2));
%     plot(regressor*0.5,'k');
% end
%
% s=6;
% hs(3)=subplot(5,2,4); hold on
% for e=1:3
%     sigvox=find(issig{e,s});
%     tcs=allts{e,s}(:,sigvox);
%     plot(mean(tcs,2));
%     plot(regressor*0.5,'k');
% end
%
% print -dpng ../figures/tmp

%%
% check echo dependence of betas
%
% s=6;
% b1=allBetas{1,s}(issig{1,s});
% b2=allBetas{2,s}(issig{1,s});
% b3=allBetas{3,s}(issig{1,s});
% figure;
% hold on
% plot(ones(1,numel(b1)),b1,'o');
% plot(2*ones(1,numel(b2)),b2,'o');
%plot(3*ones(1,numel(b3)),b3,'o');

