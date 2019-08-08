% 05/30/18
% compare grand mean BOLD in true versus control ROI
clear all; close all; clc
TR=2.8; nTR=645;

cdataFilename='../data/precomputed/allControlBoldsRoi-Filter1-Derivative1-WhiteMatter1-CSF0-LP0-30-May-2018.mat';
load(cdataFilename,'alloBolds');
allcBolds=alloBolds;

dataFilename='../data/precomputed/allBoldsRoi-Filter1-Derivative1-WhiteMatter1-CSF0-LP0-29-May-2018.mat';
load(dataFilename,'alloBolds');

allcBolds2D=cat(2,allcBolds{:});
allBolds2D=cat(2,alloBolds{:});
figure;
subplot(211); hold on
time=(0:nTR-1)*TR;
plot(time,mean(allBolds2D,2));
plot(time,mean(allcBolds2D,2));
hlg=legend('ROI','Control'); set(hlg,'box','off'); set(hlg,'orientation','horizontal');
set(hlg,'location','northwest');
print -dpng ../figures/roiVscontrolComparison
% %%
% % look at each subject's roi
% figure
% for s=1:numel(alloBolds)
%     tBold(:,s)=mean(alloBolds{s},2);
%     plot(tBold(:,s),'k'); title(s);
%     pause
% end
% 
% 
% %%
% % dataFilename='../data/precomputed/allBoldsRoi-Filter1-Derivative1-WhiteMatter1-CSF1-18-May-2018.mat';
% % load(dataFilename,'alloBolds');
% % 
% % allBolds2D=cat(2,alloBolds{:});
% % subplot(212);
% % plot(time,mean(allBolds2D,2),'k');
