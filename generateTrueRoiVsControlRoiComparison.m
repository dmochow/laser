% TMII Symposium 04/05/18
clear all; close all; clc
addpath(genpath('~/PROJECTS/COMMON'));

precomputedFilename=['../data/precomputed/allBoldsRoi ' date];
load(precomputedFilename,'allBoldsRoi','allOnsets','TR','TEs','subjStr','nSubjects','nEchos');

%%
nPerms=10;
for s=1:nSubjects
   for e=1:nEchos
        thisBold=squeeze(allBoldsRoi{s}(e,:,:));
        for p=1:nPerms
            thisBoldPerm=surrogateResponseGenerator(thisBold);
        end
        
   end
end


%%
% statistical testing
allBolds=cat(3,allBoldsRoi{:});
allControlBolds=cat(3,allBoldsControlRoi{:});
nTRs=size(allBolds,2);
nRoiVoxels=size(allBolds,3);
nControlRoiVoxels=size(allControlBolds,3);


% % normalize bolds to 0 at start time of laser
% muPre=mean(allBolds(:,1:startTR,:),2);
% allBoldsNorm=allBolds-repmat(muPre,[1 nTRs 1]);
% 
% muPreControl=mean(allControlBolds(:,1:startTR,:),2);
% allControlBoldsNorm=allControlBolds-repmat(muPreControl,[1 nTRs 1]);
% 
% [sems,mus,stds]=nansem(allBoldsNorm,3);
% [semsc,musc,stdsc]=nansem(allControlBoldsNorm,3);

%%
% fuck it, try a spectrogram
clear allS allControlS
winlen=16;
noverlap=8;
nfft=64;
fs=1/TR;
for e=1:nEchos
    for v=1:nRoiVoxels
        [s,w,t] = spectrogram(squeeze(allBolds(e,:,v)),hamming(winlen),noverlap);
        allS(:,:,e,v)=s;
    end
end

for e=1:nEchos
    for v=1:nControlRoiVoxels
        [s,w,t] = spectrogram(squeeze(allControlBolds(e,:,v)),hamming(winlen),noverlap);
        allControlS(:,:,e,v)=s;
    end
end
%%

% incoherent average
muAllS=mean(abs(allS),4);
muAllControlS=mean(abs(allControlS),4);

% coherent average
% muAllS=abs(mean(allS,4));
% muAllControlS=abs(mean(allControlS,4));

%%
figure;
for e=1:nEchos
    subplot(3,2,e*2-1);
    imagesc(muAllS(:,:,e)); cl=caxis;
    subplot(3,2,e*2);
    imagesc(muAllControlS(:,:,e)); caxis(cl);
end

figure;
for e=1:nEchos
    subplot(3,1,e);
    imagesc((muAllS(:,:,e)-muAllControlS(:,:,e))./muAllControlS(:,:,e)); colorbar
end

% %%
% figure; hold on
% plot(squeeze(mean(muAllS(1:10,:,3),1)));
% plot(squeeze(mean(muAllControlS(1:10,:,3),1)));
% %%
% nTRs=size(allBolds,2);
% pval=zeros(nEchos,nTRs);
% alpha=0.01;
% for e=1:nEchos
%     for t=1:nTRs
%         pval(e,t)=ranksum(squeeze(allBoldsNorm(e,t,:)),squeeze(allControlBoldsNorm(e,t,:)),'tail','right');
%     end
% end
% [pfdr,isSig]=fdr(pval,alpha);
% 
% % figure
% time=(0:644)*TR;
% tLaser=[600 1200];
% 
% figure;
% for e=1:nEchos
%     subplot(nEchos,1,e); hold on
%     shadedErrorBar(time,mus(e,:),sems(e,:),'b',1); title('Echo 1');
%     shadedErrorBar(time,musc(e,:),semsc(e,:),'r',1); title('Echo 1');     yl=ylim;
%     plot(time(isSig(e,:)),yl(1)*ones(sum(isSig(e,:)),1),'*k','LineStyle','none');
% 
%     harea=area(tLaser,yl(2)*[1 1],'BaseValue',yl(1));
%     set(harea,'FaceColor',[0.7 0.7 0.7]);
%     set(harea,'FaceAlpha',0.25);
%     set(harea,'EdgeColor','none');
% end
% 
% %%
% % grand mean in control roi
% 
% figure;
% subplot(411); hold on
% plot(time,grandControlBold(1,:)); title('Echo 1');
% yl=ylim;
% harea=area(tLaser,yl(2)*[1 1],'BaseValue',yl(1));
% set(harea,'FaceColor',[0.7 0.7 0.7]);
% set(harea,'FaceAlpha',0.5);
% set(harea,'EdgeColor','none');
% subplot(412); hold on
% plot(time,grandControlBold(2,:)); title('Echo 2');
% yl=ylim;
% harea=area(tLaser,yl(2)*[1 1],'BaseValue',yl(1));
% set(harea,'FaceColor',[0.7 0.7 0.7]);
% set(harea,'FaceAlpha',0.5);
% set(harea,'EdgeColor','none');
% subplot(413); hold on
% plot(time,grandControlBold(3,:)); title('Echo 3');
% yl=ylim;
% harea=area(tLaser,yl(2)*[1 1],'BaseValue',yl(1));
% set(harea,'FaceColor',[0.7 0.7 0.7]);
% set(harea,'FaceAlpha',0.5);
% set(harea,'EdgeColor','none');
% subplot(414); hold on
% plot(time,mean(grandControlBold)); title('Echo mean');
% yl=ylim;
% harea=area(tLaser,yl(2)*[1 1],'BaseValue',yl(1));
% set(harea,'FaceColor',[0.7 0.7 0.7]);
% set(harea,'FaceAlpha',0.5);
% set(harea,'EdgeColor','none');
% print('-dpng','-r600',['../figures/grandMeanControltimeCourse ' date]);
% 
% %%
% 
% %% 
% % % show subject means
% % sBold=zeros(nEchos,nTRs,nSubjects);
% % for s=1:nSubjects, sBold(:,:,s)=mean(allBoldsRoi{s},3); end
% % musBold=mean(sBold,3);
% % figure
% % for e=1:nEchos
% %     subplot(nEchos,1,e);
% %     plot(squeeze(sBold(e,:,:)));
% %     hold on
% %     plot(musBold(e,:),'k','LineWidth',4);
% % end
% % 
% % %%
% % figure
% % for s=1:nSubjects
% %     plot(sBold(2,:,s),'k');
% %     title(subjStr{s});
% %     pause
% % end
