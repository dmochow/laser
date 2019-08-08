% 06/04/18
% bench for testing chow test
clear all; close all; clc
TR=2.8; nTR=645; time=(0:nTR-1)*TR; 
laserOnsetTR=215;
laserOffsetTR=430; 
dataFilename='../data/precomputed/allBoldsRoi-Filter1-Derivative1-WhiteMatter1-CSF0-LP0-STANDARDIZE1-TSHIFT1-Echo1-05-Jun-2018.mat';
%dataFilename='../data/precomputed/allBoldsRoi-Filter1-Derivative1-WhiteMatter1-CSF0-LP0-STANDARDIZE1-Echo2-31-May-2018.mat';
%dataFilename='../data/precomputed/allBoldsRoi-Filter1-Derivative1-WhiteMatter1-CSF0-LP0-STANDARDIZE1-Echo3-31-May-2018.mat';

load(dataFilename,'alloBolds');
allBolds2D=cat(2,alloBolds{:});
[semBold,grandBold]=nansem(allBolds2D,2);

nlags=5; % optimized in optimizeARmodelOrder.m

%% Analysis 1: compare pre with peri, pre with post
y=grandBold;
X=tplitz(y,nlags); X=X(:,2:end-1);

% is pre different from peri?
[h1,p1] = chowtest(X(1:laserOffsetTR,:),y(1:laserOffsetTR),laserOnsetTR);

% is pre different from post?
[h2,p2] = chowtest(cat(1,X(1:laserOnsetTR,:),X(laserOffsetTR+1:end,:)),cat(1,y(1:laserOnsetTR),y(laserOffsetTR+1:end)),laserOnsetTR);

% [pfdr,isSig]=fdr(pvals,alpha);
% notSig=~isSig;
% %% draw figure
% figure;
% hs(1)=subplot(211);
% hsh=shadedErrorBar(time,grandBold,semBold);  
% hs(2)=subplot(212); hold on
% hp=plot(bpsTime,pvals);
% %plot(bpsTime(isSig),pvals(isSig),'*b','MarkerFaceColor','b');
% plot(bpsTime(notSig),pvals(notSig),'or','MarkerFaceColor','r');
% ylabel('P-value');
% xlabel('Time (s)');
% 
% 
% %
% set(get(hs(1),'ylabel'),'String','BOLD (a.u.)');
% set(hs(1),'Xtick',[0 600 1200 1800]);
% set(hs(1),'XTickLabel',{'0:00','10:00','20:00','30:00'});
% set(hs(2),'Xtick',[0 600 1200 1800]);
% set(hs(2),'XTickLabel',{'0:00','10:00','20:00','30:00'});

%% Analysis 3: split at 2-minute windows following laser onset
%% set up model and data structures
nlags=5;
bps2=laserOnsetTR;
bpsTime2=bps2*TR;
wins=[(10:1:29)' (11:1:30)'];
wins=round(wins*60/TR);
winsTime=wins*TR;
nWins=size(wins,1);
pvals2=zeros(nWins,1);
for w=1:nWins
    w
    y=cat(1,grandBold(1:laserOnsetTR),grandBold(wins(w,1):wins(w,2)));
    %y=surrogateResponseGenerator(y);
    X=tplitz(y,nlags); X=X(:,2:end-1);
    [h,p] = chowtest(X,y,bps2);
    pvals2(w)=p;
end
% run Chow tests

alpha=0.05;
isSig2=pvals2<alpha;
[pfdr2,isSig2]=fdr(pvals2,alpha);
% notSig=~isSig;
% 
% % draw figure
% figure;
% hs(1)=subplot(211);
% hsh=shadedErrorBar(time,grandBold,semBold);  
% hs(2)=subplot(212); 
% hp=plot(winsTime(:,1),pvals2,'*k'); hold on
% plot(winsTime(isSig2),pvals2(isSig2),'*b','MarkerFaceColor','b');
% %plot(bpsTime(notSig),pvals(notSig),'or','MarkerFaceColor','r');
% ylabel('P-value');
% xlabel('Time (s)');
% set(hs(1),'Xtick',[0 600 1200 1800]);
% set(hs(2),'Xtick',[0 600 1200 1800]);
% set(hs(2),'Xlim',[0 1800]); set(hs(2),'Ylim',[0 1]);
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
