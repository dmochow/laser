% compare effect size across echos, evaluate whether it's linear
clear all; close all; clc

%%
PRINT_FIGURE=1;
pathToData='../data/SAVG/NII';
anatFilename=fullfile(pathToData,'TT_N27+tlrc');
brainMaskFilename=fullfile(pathToData,'resampled_brain_mask+tlrc.BRIK');
roiMaskFilename=fullfile(pathToData,'mu_roi_r19_z26+tlrc.BRIK');
[~, brainMask, ~, ~] = BrikLoad (brainMaskFilename);
brainMask=logical(brainMask);

colors=[0 0.4470 0.7410;0.8500    0.3250    0.0980];
for e=1:3
    isSigBefDurFilename=fullfile(pathToData,['thresh_isSigBefDur_NLAGS5_CONSTANT0_s8muwBoldEcho' num2str(e) '-Filter1-Derivative1-WhiteMatter3-STANDARDIZE1-TSHIFT1-CSHIFT2+tlrc.BRIK']);
    tstatsBefDurFilename=fullfile(pathToData,['tstatsBefDur_NLAGS5_CONSTANT0_s8muwBoldEcho' num2str(e) '-Filter1-Derivative1-WhiteMatter3-STANDARDIZE1-TSHIFT1-CSHIFT2+tlrc.BRIK']);
    
    
    %% load the anatomy and the significance clusters
    [~, isSig, ~, ~] = BrikLoad (isSigBefDurFilename);
    [~, tstats, ~, ~] = BrikLoad (tstatsBefDurFilename);
    
    allIsSig(:,e)=vol2ts(logical(isSig),brainMask);
    alltstats(:,e)=vol2ts(tstats,brainMask);
    
end

%%
xaxis=[12.8,34.3,55.6]';
figure;
hs(1)=subplot(221);
%hp(:,1)=plot(xaxis,alltstats(allIsSig(:,1),:).','o-','MarkerFaceColor','auto');
hp(:,1)=plot(xaxis,alltstats(allIsSig(:,1),:).','o-');
hs(2)=subplot(222);
hp(:,2)=plot(xaxis,mean(alltstats(allIsSig(:,1),:).',2),'o-');
%hp(:,2)=plot(xaxis,mean(alltstats(allIsSig(:,1),:).',2),'o-','MarkerFaceColor','auto');
for l=1:size(hp,1)
    clr=get(hp(l,1),'Color');
    set(hp(l,1),'MarkerFaceColor',clr);
    
    clr=get(hp(l,2),'Color');
    set(hp(l,2),'MarkerFaceColor',clr);
end
%%
% esthetics
set(hs(:),'XLim',[xaxis(1)-3 xaxis(end)+3]);
set(hs(:),'YLim',[0 10]);
set(hs(:),'Xtick',xaxis);
set(hs(:),'Ytick',0:2:10);
set(get(hs(1),'YLabel'),'String','F statistic','FontSize',16)
set(get(hs(1),'XLabel'),'String','Echo time (ms)','FontSize',16)
set(get(hs(2),'XLabel'),'String','Echo time (ms)','FontSize',16)
%set(get(hs(1),'Title'),'String','All significant voxels','FontSize',16,'FontWeight','normal')
%set(get(hs(2),'Title'),'String','Mean','FontSize',16,'FontWeight','normal')
set(hs(:),'Box','off');
%%
% sublabel
if PRINT_FIGURE
    sublabel([hs(1) hs(2)],-20,-20,'FontSize',16,'FontWeight','Bold');
    print -dpng ../figures/effectVsEcho
    crop('../figures/effectVsEcho.png');
end

%%
% stats
x=xaxis;
y=mean(alltstats(allIsSig(:,1),:).',2);
Xr=[ones(3,1)];
Xf=[x ones(3,1)];
Jf=(eye(3)-Xf*pinv(Xf))*y;
Jr=(eye(3)-Xr*pinv(Xr))*y;
SSf=sum(Jf.^2);
SSr=sum(Jr.^2 );
DFf=3-2;
DFr=3-1;
F= ( (SSr-SSf)/(DFr-DFf) ) / (SSf/DFf)
p = 1-fcdf(F,DFr,DFf)
r2=1-SSf/SSr

% y=Bo vs y = B*x + Bo