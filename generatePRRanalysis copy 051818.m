% TMII Symposium 04/05/18
% BRASH 04/14/18
% TMII poster 04/22/18
% paper 05/18/18
clear all; close all; clc
addpath(genpath('~/PROJECTS/COMMON'));

%precomputedFilename=['../data/precomputed/allBoldsRoi ' date];
%load(precomputedFilename,'allBoldsRoiOut','allOnsets','TR','TEs','subjStr','nSubjects','nEchos');
precomputedFilename='../data/precomputed/allBoldsRoi-Filter1-Derivative1-WhiteMatter1-16-May-2018.mat';
load(precomputedFilename,'alloBolds');
nEchos=1; nSubjects=numel(alloBolds);
%%
nTRs=645;
nPerms=500;
%muBoldsPerm=zeros(nEchos,nTRs,nPerms);
muBoldsPerm=zeros(nTRs,nPerms);
for p=1:nPerms
    p
    allBoldsRoiPerm=cell(nSubjects,1);
    for s=1:nSubjects
        for e=1:nEchos
            %thisBold=squeeze(allBoldsRoiOut{s}(e,:,:));
            %allBoldsRoiPerm{s}(e,:,:)=surrogateResponseGenerator(thisBold);
            thisBold=squeeze(alloBolds{s}(:,:));
            allBoldsRoiPerm{s}=surrogateResponseGenerator(thisBold);
        end  
        tmp=cat(3,allBoldsRoiPerm{:});
        muBoldsPerm(:,:,p)=mean(tmp,3);
    end
end


%%
muBolds=mean(cat(3,allBoldsRoiOut{:}),3);
pval=mean(muBoldsPerm>repmat(muBolds,[1 1 nPerms]),3);
pval(:,1:214)=NaN;
%pval(:,429:end)=NaN;
alpha=0.05;
[pfdr,isSig]=fdr(pval,alpha);
isSig=pval<alpha;

%%
% only keep time points that at least appear consecutively
isSig2=zeros(size(isSig));
for e=1:nEchos
    tmp=find(isSig(e,:));
    p=find(diff(tmp)==1);
    q=[p;p+1];
    tmp2=tmp(q);
    isSig2(e,unique(tmp2(:)))=1;
end
%%

time=(0:nTRs-1)*TR;
tLaser=[600 1200];
figure;
for e=1:nEchos
    hs(e)=subplot(nEchos,1,e); hold on
    plot(time,muBolds(e,:),'b'); 
    title(['Echo ' num2str(e)]);
    ylim([-50 50]);
    yl=ylim;
    hsig=plot(time(isSig2(e,:)==1),yl(1)*ones(sum(isSig2(e,:)==1),1),'*k','LineStyle','none');

    harea=area(tLaser,yl(2)*[1 1],'BaseValue',yl(1));
    set(harea,'FaceColor',[0.7 0.7 0.7]);
    set(harea,'FaceAlpha',0.25);
    set(harea,'EdgeColor','none');
    
    if e==1
        hlg=legend([harea hsig],'Laser On','p<0.05');
    end
end
set(get(hs(3),'xlabel'),'string','Time (s)');
set(get(hs(2),'ylabel'),'string','BOLD (a.u.)');
set(hlg,'box','off');
set(hlg,'orientation','horizontal');
print('-dpng',['../figures/PRR_analysis' date]);

