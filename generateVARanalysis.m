% 04/23/18
% compare the BOLD in the ROI with mock BOLD generated by forecasting from
% the first 10 minutes using vector autoregressive model forecasting
clear all; close all; clc
addpath(genpath('~/PROJECTS/COMMON'));

%%
precomputedFilename='../data/precomputed/allBoldsRoi 22-Apr-2018';
load(precomputedFilename,'allBoldsRoiOut','allOnsets','TR','TEs','subjStr','nSubjects','nEchos');

%%
nTRs=645;
nPerms=5;
p=3; % VAR model order
nVoxelsMax=7;
allBoldsRoiPerm=cell(nSubjects,1);
vxlRange=10:10+nVoxelsMax-1;
for s=1:nSubjects
    for e=1:nEchos
        [s,e]
        thisBold=squeeze(allBoldsRoiOut{s}(e,:,vxlRange));
        thisBoldPre=thisBold(1:allOnsets(s,1),:); % TODO: use actual onsets
        try
            thisBoldForecast = fitAndForecast(thisBoldPre,p,nTRs-allOnsets(s,1),nPerms);
            thisBoldForecast=cat(1,repmat(thisBoldPre,[1 1 nPerms]),thisBoldForecast);
            allBoldsRoiPerm{s}(e,:,:,:)=thisBoldForecast;
        catch
            warning('Something went wrong'); 
            allBoldsRoiPerm{s}(e,:,:,:)=zeros(nTRs,nVoxelsMax,nPerms);
        end
    end
end
tmp=cat(3,allBoldsRoiPerm{:});
muBoldsPerm=squeeze(mean(tmp,3));


% %%
% % real bold
% for s=1:nSubjects
%     ctmp{s}=allBoldsRoiOut{s}(:,:,vxlRange);
% end
% muBolds=mean(cat(3,ctmp{:}),3);
% 
% %%
% % plot
% figure; 
% hold on
% plot(muBolds.','k','LineWidth',4);
% plot(muBoldsPerm(:,:,1).')




% pval=mean(muBoldsPerm>repmat(muBolds,[1 1 nPerms]),3);
% pval(:,1:214)=NaN;
% %pval(:,429:end)=NaN;
% alpha=0.05;
% [pfdr,isSig]=fdr(pval,alpha);
% isSig=pval<alpha;
% 
% %%
% % only keep time points that at least appear consecutively
% isSig2=zeros(size(isSig));
% for e=1:nEchos
%     tmp=find(isSig(e,:));
%     p=find(diff(tmp)==1);
%     q=[p;p+1];
%     tmp2=tmp(q);
%     isSig2(e,unique(tmp2(:)))=1;
% end
% %%
% 
% time=(0:nTRs-1)*TR;
% tLaser=[600 1200];
% figure;
% for e=1:nEchos
%     hs(e)=subplot(nEchos,1,e); hold on
%     plot(time,muBolds(e,:),'b'); 
%     title(['Echo ' num2str(e)]);
%     ylim([-50 50]);
%     yl=ylim;
%     hsig=plot(time(isSig2(e,:)==1),yl(1)*ones(sum(isSig2(e,:)==1),1),'*k','LineStyle','none');
% 
%     harea=area(tLaser,yl(2)*[1 1],'BaseValue',yl(1));
%     set(harea,'FaceColor',[0.7 0.7 0.7]);
%     set(harea,'FaceAlpha',0.25);
%     set(harea,'EdgeColor','none');
%     
%     if e==1
%         hlg=legend([harea hsig],'Laser On','p<0.05');
%     end
% end
% set(get(hs(3),'xlabel'),'string','Time (s)');
% set(get(hs(2),'ylabel'),'string','BOLD (a.u.)');
% set(hlg,'box','off');
% set(hlg,'orientation','horizontal');
% print('-dpng',['../figures/PRR_analysis' date]);
% 
