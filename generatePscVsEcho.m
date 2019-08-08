% 10.29.18
% running this to get PSC and regress onto echo
% do it on all 214 significant voxels
clear all; close all; clc
PRINT_FIGURE=1;
nTRs=645;
TR=2.8;
fs=1/TR; % for AR spectra
nfft=32; % for AR spectra
freqs=(0:nfft-1)/nfft*fs;
posFreqMask=1:nfft/2+1;
xaxis=[12.8,34.3,55.6]';
colors=[0 0.4470 0.7410;0.8500    0.3250    0.0980];
pathToData='../data/SAVG/NII';
anatFilename=fullfile(pathToData,'TT_N27+tlrc');
muBoldFilenames={'s8muwBoldEcho1-Filter1-Derivative1-WhiteMatter3-STANDARDIZE1-TSHIFT1-CSHIFT2+tlrc.BRIK',...
    's8muwBoldEcho2-Filter1-Derivative1-WhiteMatter3-STANDARDIZE1-TSHIFT1-CSHIFT2+tlrc.BRIK',...
    's8muwBoldEcho3-Filter1-Derivative1-WhiteMatter3-STANDARDIZE1-TSHIFT1-CSHIFT2+tlrc.BRIK'};
filenames={'thresh_isSigBefDur_NLAGS5_s8muwBoldEcho1-Filter1-Derivative1-WhiteMatter3-STANDARDIZE1-TSHIFT1-CSHIFT2+tlrc.BRIK',...
    'thresh_isSigBefDur_NLAGS5_s8muwBoldEcho2-Filter1-Derivative1-WhiteMatter3-STANDARDIZE1-TSHIFT1-CSHIFT2+tlrc.BRIK',...
    'thresh_isSigBefDur_NLAGS5_s8muwBoldEcho3-Filter1-Derivative1-WhiteMatter3-STANDARDIZE1-TSHIFT1-CSHIFT2+tlrc.BRIK'};
arFilenames={'arCoeffsResidsEcho1-NLAGS5_CONSTANT0.mat','arCoeffsResidsEcho2-NLAGS5_CONSTANT0.mat','arCoeffsResidsEcho3.mat'};
nEchos=3;

brainMaskFilename='resampled_brain_mask+tlrc.BRIK';
[~,brainMask,~,~]=BrikLoad(fullfile(pathToData,brainMaskFilename));
brainMask=logical(brainMask);

[XX,YY,ZZ]=ndgrid(1:size(brainMask,1),1:size(brainMask,2),1:size(brainMask,3));

hf=figure;
for e=1:nEchos
    
    % get significant voxels
    thisBefDurFilename=['thresh_isSigBefDur_s8muwBoldEcho1-Filter1-Derivative1-WhiteMatter3-STANDARDIZE1-TSHIFT1-CSHIFT2+tlrc.BRIK'];
    [~,isSigBefDur1,~,~]=BrikLoad(fullfile(pathToData,thisBefDurFilename));
    
    thisBefDurFilename=['thresh_isSigBefDur_s8muwBoldEcho2-Filter1-Derivative1-WhiteMatter3-STANDARDIZE1-TSHIFT1-CSHIFT2+tlrc.BRIK'];
    [~,isSigBefDur2,~,~]=BrikLoad(fullfile(pathToData,thisBefDurFilename));
    
    thisBefDurFilename=['thresh_isSigBefDur_s8muwBoldEcho3-Filter1-Derivative1-WhiteMatter3-STANDARDIZE1-TSHIFT1-CSHIFT2+tlrc.BRIK'];
    [~,isSigBefDur3,~,~]=BrikLoad(fullfile(pathToData,thisBefDurFilename));
    
    
    %isSig1DbefDur=vol2ts(isSigBefDur,brainMask);
    isSigBefDur=(isSigBefDur1|isSigBefDur2|isSigBefDur3) & YY<32;
    isSig1DbefDur=vol2ts(isSigBefDur,brainMask);
    
    % AR coefficients
    load(fullfile(pathToData,arFilenames{e}),'coeffsPreStim','coeffsPrePost');
    
    % before vs during
    arMagSpec = getARspectrum(coeffsPreStim,nfft);
    arMagSpecSig=arMagSpec(:,isSig1DbefDur>0,:);
    
    % 10.16.18 percent signal change
    tmp1=squeeze(arMagSpecSig(posFreqMask,:,1));
    tmp2=squeeze(sum(arMagSpecSig(posFreqMask,:,1)));
    freqWeighting=tmp1./repmat(tmp2,size(tmp1,1),1);
    
    pscf=( arMagSpecSig(posFreqMask,:,2)-arMagSpecSig(posFreqMask,:,1) )./arMagSpecSig(posFreqMask,:,1);
    psc(:,e)=sum(freqWeighting .* pscf,1);
   
end



x=xaxis;
y=mean(psc,1)';
Xr=ones(3,1);
Xf=[x ones(3,1)];
Jf=(eye(3)-Xf*pinv(Xf))*y;
Jr=(eye(3)-Xr*pinv(Xr))*y;
SSf=sum(Jf.^2);
SSr=sum(Jr.^2 );
DFf=3-2;
DFr=3-1;
F= ( (SSr-SSf)/(DFr-DFf) ) / (SSf/DFf);
pval = 1-fcdf(F,DFr,DFf);
r2=1-SSf/SSr;


%%
[B,BINT] = regress(mean(psc,1)'*100,[xaxis ones(3,1)]) ;


figure;
hs(1)=subplot(221);
hp(:,1)=plot(xaxis,psc*100,'o-','MarkerSize',6); yl=ylim;
hs(2)=subplot(222);
hp(:,2)=plot(xaxis,mean(psc,1)*100,'o','LineStyle','none','MarkerSize',6); hold on
%ylim(yl);
xl=xlim;
hreg=plot(linspace(xl(1),xl(2),100),B(2)+B(1)*linspace(xl(1),xl(2),100));
set(hreg,'Color','b');
%set(hreg)
%plot
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
%set(hs(:),'YLim',[0 10]);
set(hs(:),'Xtick',xaxis);
%set(hs(:),'Ytick',0:2:10);
set(get(hs(1),'YLabel'),'String','Percent Change','FontSize',16)
set(get(hs(1),'XLabel'),'String','Echo time (ms)','FontSize',16)
set(get(hs(2),'XLabel'),'String','Echo time (ms)','FontSize',16)
%set(get(hs(1),'Title'),'String','All significant voxels','FontSize',16,'FontWeight','normal')
%set(get(hs(2),'Title'),'String','Mean','FontSize',16,'FontWeight','normal')
set(hs(:),'Box','off');

%%
pos1=get(hs(1),'Position');
set(hs(1),'Position',[0.071 pos1(2) pos1(3) pos1(4)]);

pos2=get(hs(2),'Position');
set(hs(2),'Position',[0.44 pos2(2) pos2(3) pos2(4)]);
%%
% sublabel
if PRINT_FIGURE
    sublabel([hs(1) hs(2)],-15,-30,'FontSize',16,'FontWeight','Bold');
    print -dpng ../figures/pscVsEcho
    crop('../figures/pscVsEcho.png');
end

%%
% % compare effect size across echos, evaluate whether it's linear
% clear all; close all; clc
%
% %%
% PRINT_FIGURE=1;
% pathToData='../data/SAVG/NII';
% anatFilename=fullfile(pathToData,'TT_N27+tlrc');
% brainMaskFilename=fullfile(pathToData,'resampled_brain_mask+tlrc.BRIK');
% roiMaskFilename=fullfile(pathToData,'mu_roi_r19_z26+tlrc.BRIK');
% [~, brainMask, ~, ~] = BrikLoad (brainMaskFilename);
% brainMask=logical(brainMask);
%
% colors=[0 0.4470 0.7410;0.8500    0.3250    0.0980];
% for e=1:3
%     isSigBefDurFilename=fullfile(pathToData,['thresh_isSigBefDur_NLAGS5_CONSTANT0_s8muwBoldEcho' num2str(e) '-Filter1-Derivative1-WhiteMatter3-STANDARDIZE1-TSHIFT1-CSHIFT2+tlrc.BRIK']);
%     tstatsBefDurFilename=fullfile(pathToData,['tstatsBefDur_NLAGS5_CONSTANT0_s8muwBoldEcho' num2str(e) '-Filter1-Derivative1-WhiteMatter3-STANDARDIZE1-TSHIFT1-CSHIFT2+tlrc.BRIK']);
%
%
%     %% load the anatomy and the significance clusters
%     [~, isSig, ~, ~] = BrikLoad (isSigBefDurFilename);
%     [~, tstats, ~, ~] = BrikLoad (tstatsBefDurFilename);
%
%     allIsSig(:,e)=vol2ts(logical(isSig),brainMask);
%     alltstats(:,e)=vol2ts(tstats,brainMask);
%
% end
%
% %%
% xaxis=[12.8,34.3,55.6]';
% figure;
% hs(1)=subplot(221);
% %hp(:,1)=plot(xaxis,alltstats(allIsSig(:,1),:).','o-','MarkerFaceColor','auto');
% hp(:,1)=plot(xaxis,alltstats(allIsSig(:,1),:).','o-');
% hs(2)=subplot(222);
% hp(:,2)=plot(xaxis,mean(alltstats(allIsSig(:,1),:).',2),'o-');
% %hp(:,2)=plot(xaxis,mean(alltstats(allIsSig(:,1),:).',2),'o-','MarkerFaceColor','auto');
% for l=1:size(hp,1)
%     clr=get(hp(l,1),'Color');
%     set(hp(l,1),'MarkerFaceColor',clr);
%
%     clr=get(hp(l,2),'Color');
%     set(hp(l,2),'MarkerFaceColor',clr);
% end
% %%
% % esthetics
% set(hs(:),'XLim',[xaxis(1)-3 xaxis(end)+3]);
% set(hs(:),'YLim',[0 10]);
% set(hs(:),'Xtick',xaxis);
% set(hs(:),'Ytick',0:2:10);
% set(get(hs(1),'YLabel'),'String','F statistic','FontSize',16)
% set(get(hs(1),'XLabel'),'String','Echo time (ms)','FontSize',16)
% set(get(hs(2),'XLabel'),'String','Echo time (ms)','FontSize',16)
% %set(get(hs(1),'Title'),'String','All significant voxels','FontSize',16,'FontWeight','normal')
% %set(get(hs(2),'Title'),'String','Mean','FontSize',16,'FontWeight','normal')
% set(hs(:),'Box','off');
% %%
% % sublabel
% if PRINT_FIGURE
%     sublabel([hs(1) hs(2)],-20,-20,'FontSize',16,'FontWeight','Bold');
%     print -dpng ../figures/effectVsEcho
%     crop('../figures/effectVsEcho.png');
% end
%
% %%
% % stats
% x=xaxis;
% y=mean(alltstats(allIsSig(:,1),:).',2);
% Xr=[ones(3,1)];
% Xf=[x ones(3,1)];
% Jf=(eye(3)-Xf*pinv(Xf))*y;
% Jr=(eye(3)-Xr*pinv(Xr))*y;
% SSf=sum(Jf.^2);
% SSr=sum(Jr.^2 );
% DFf=3-2;
% DFr=3-1;
% F= ( (SSr-SSf)/(DFr-DFf) ) / (SSf/DFf)
% p = 1-fcdf(F,DFr,DFf)
% r2=1-SSf/SSr
%
% % y=Bo vs y = B*x + Bo