% TMII Symposium 04/05/18
% BRASH 04/14/18
% TMII poster 04/22/18
% paper 05/18/18
% paper 06/22/18
clear all; close all; clc
addpath(genpath('~/PROJECTS/COMMON'));

precomputedFilename='../data/precomputed/allBoldsRoi-Filter1-Derivative1-WhiteMatter1-CSF0-LP0-STANDARDIZE1-TSHIFT1-Echo1-05-Jun-2018';
load(precomputedFilename,'alloBolds');
nSubjects=numel(alloBolds);
TR=2.8;
nTRs=645;
nPerms=1000;
firstTRkeep=3; % ignore first few unstable frames
muBoldsPerm=zeros(nTRs-(firstTRkeep-1),nPerms);
for p=1:nPerms
    p
    allBoldsRoiPerm=cell(nSubjects,1);
    for s=1:nSubjects
        thisBold=squeeze(alloBolds{s}(:,:));
        thisBold=thisBold(firstTRkeep:end,:); % ignore first few TRs
        allBoldsRoiPerm{s}=surrogateResponseGenerator(thisBold);
        tmp=cat(2,allBoldsRoiPerm{:});
        muBoldsPerm(:,p)=mean(tmp,2);
    end
end


%%
muBolds=mean(cat(2,alloBolds{:}),2); % the true grand-mean
muBolds=muBolds(firstTRkeep:end);
pval=mean(muBoldsPerm>repmat(muBolds,[1 nPerms]),2);
%pval(:,1:214)=NaN;
%pval(:,429:end)=NaN;
alpha=0.05;
[pfdr,isSig]=fdr(pval,alpha);
isSig=pval<alpha;

%%
% only keep time points that at least appear consecutively
isSig2=zeros(size(isSig));
%for e=1:nEchos
    tmp=find(isSig);
    p=find(diff(tmp)==1);
    q=[p;p+1];
    tmp2=tmp(q);
    isSig2(unique(tmp2(:)))=1;
%end
%%

time=(firstTRkeep-1:nTRs-1)*TR;
tLaser=[600 1200];
figure;
%for e=1:nEchos
    hs(1)=subplot(2,1,1); hold on
    plot(time,muBolds(),'b');
    title(['Echo ' num2str(1)]);
    %ylim([-50 50]);
    yl=ylim;
    hsig=plot(time(isSig2==1),yl(1)*ones(sum(isSig2==1),1),'*k','LineStyle','none');
    
    harea=area(tLaser,yl(2)*[1 1],'BaseValue',yl(1));
    set(harea,'FaceColor',[0.7 0.7 0.7]);
    set(harea,'FaceAlpha',0.25);
    set(harea,'EdgeColor','none');
    
    %if e==1
        hlg=legend([harea hsig],'Laser On','p<0.05');
    %end
%end
set(get(hs(1),'xlabel'),'string','Time (s)');
set(get(hs(1),'ylabel'),'string','BOLD (a.u.)');
set(hlg,'box','off');
set(hlg,'orientation','horizontal');
print('-dpng',['../figures/PRR_analysis' date]);

