% 06/04/18
% bench for testing chow test
clear all; close all; clc
TR=2.8; nTR=645; time=(0:nTR-1)*TR; 
laserOnsetTR=215;
laserOffsetTR=430; 
nlags=5; % optimized in optimizeARmodelOrder.m
alpha=0.05;
dataFilename='../data/precomputed/allBoldsRoi-Filter1-Derivative1-WhiteMatter1-CSF0-LP0-STANDARDIZE1-TSHIFT1-Echo1-05-Jun-2018.mat';
dataFilename='../data/precomputed/allControlBoldsRoi-Filter1-Derivative1-WhiteMatter1-CSF0-LP0-STANDARDIZE1-TSHIFT1-Echo1-04-Jun-2018';
load(dataFilename,'alloBolds');
allBolds2D=cat(2,alloBolds{:});
[semBold,grandBold]=nansem(allBolds2D,2);

bps2=laserOnsetTR;
bpsTime2=bps2*TR;
winlen=2; winshift=0.25;
wins=(10:winshift:27)';
wins(:,2)=wins(:,1)+winlen;
wins=round(wins*60/TR);
winsTime=wins*TR;
nWins=size(wins,1);
pvals2=zeros(nWins,1);
for w=1:nWins
    w
    y=cat(1,grandBold(1:laserOnsetTR),grandBold(wins(w,1):wins(w,2)));
    X=tplitz(y,nlags); X=X(:,2:end-1);
    [h,p] = chowtest(X,y,bps2);
    pvals2(w)=p;
end



%isSig2=pvals2<alpha;
[pfdr2,isSig2]=fdr(pvals2,alpha);
% notSig=~isSig;
% 
% draw figure
figure;
hs(1)=subplot(211);
hsh=shadedErrorBar(time,grandBold,semBold);  
hs(2)=subplot(212); 
hp=plot(winsTime(:,1),pvals2,'*k'); hold on
plot(winsTime(isSig2),pvals2(isSig2),'*b','MarkerFaceColor','b');
%plot(bpsTime(notSig),pvals(notSig),'or','MarkerFaceColor','r');
ylabel('P-value');
xlabel('Time (s)');
set(hs(1),'Xtick',[0 600 1200 1800]);
set(hs(2),'Xtick',[0 600 1200 1800]);
set(hs(2),'Xlim',[0 1800]); set(hs(2),'Ylim',[0 1]);
% % Mdl = fitlm(X,y);
% % 
% % %
% % res = Mdl.Residuals.Raw;
% % figure;
% % subplot(221);
% % plotResiduals(Mdl,'lagged');
% % subplot(222);
% % plotResiduals(Mdl,'caseorder');
% 
