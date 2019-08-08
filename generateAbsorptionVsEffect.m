clear all; close all; clc
% attempt to correlate absorption with BOLD effect

% TODO: check to make sure that the two measures are in fact aligned, as
% one went through the Afni process (absorption) and the other didn't

%% first bring in absorption
pathToData='../data/SAVG/NII';
absorptionFilename=fullfile(pathToData,'Amu+tlrc.BRIK');
[~, absorption, ~, ~] = BrikLoad (absorptionFilename);
brainMaskFilename='brain_mask+tlrc';
[~, brainMask, info, ~] = BrikLoad (fullfile(pathToData,brainMaskFilename));
%absorption=absorption(end:-1:1,:,:);
fs=14; % font size for axis labels and text

%% in order to prune significant voxel list in echo 3
[XX,YY,ZZ]=ndgrid(1:size(brainMask,1),1:size(brainMask,2),1:size(brainMask,3));
% figure
% imagesc(brainMask(:,:,88));

% min max of sig voxels
% %1
% 59    93
% 14    31
% 74   100
% %2
% 59    86
% 16    31
% 76    95
% %3
% 64    88
% 19    98
% 86   123

% ...so take YY<32 to exluce spurious clusters
%% now bring in the volume of effect sizes (test stats and p-values)
%load('../data/SAVG/NII/tstatVolume.mat','vstat1','vstat2','vpvals1','vpvals2');

for e=1:3  % echo index
    isSigBefDurFilename=fullfile(pathToData,['resampled_thresh_isSigBefDur_s8muwBoldEcho' num2str(e) '-Filter1-Derivative1-WhiteMatter3-STANDARDIZE1-TSHIFT1-CSHIFT2+tlrc.BRIK']);
    isSigBefAftFilename=fullfile(pathToData,['resampled_thresh_isSigBefAft_s8muwBoldEcho' num2str(e) '-Filter1-Derivative1-WhiteMatter3-STANDARDIZE1-TSHIFT1-CSHIFT2+tlrc.BRIK']);
    
    [~, isSigBefDur, ~, ~] = BrikLoad (isSigBefDurFilename);
    [~, isIsSigBefAft, ~, ~] = BrikLoad (isSigBefAftFilename);
    
    % define two sets of voxels
    indsIn=brainMask(:)>0 & absorption(:)~=0 & isSigBefDur(:)==1 & YY(:)<32;
    indsOut=brainMask(:)>0 & absorption(:)~=0 & isSigBefDur(:)==0;
    % brain mask, region of support of absorption
    
    absIn=absorption(indsIn);
    absOut=absorption(indsOut);
    xIn=XX(indsIn); yIn=YY(indsIn); zIn=ZZ(indsIn);
    [min(xIn),max(xIn); min(yIn),max(yIn); min(zIn),max(zIn)]
    [e, mean(absIn), mean(absOut)]
    
    
    
    absIn=log(absIn); absOut=log(absOut);
    
    % run logistic regression and get AUC
    
    numel(absOut),numel(absIn)
    R=round(numel(absOut)/numel(absIn));
    %R=1;
    absInR=repmat(absIn,R,1);
    X=[absInR;absOut];
    Y=[1*ones(size(absInR,1),1); 2*ones(size(absOut,1),1)];
    [B,dev,stats]=mnrfit(X,Y);
    pihat = mnrval(B,X);
    
    scores=B(2)*X+B(1);
    yhat=scores>0;
    thresh=-B(1)/B(2);
    thresh95= -( log( 1/0.95-1 ) + B(1) ) / B(2);
    thresh90= -( log( 1/0.9-1 ) + B(1) ) / B(2);
    thresh50= -( 0 + B(1) ) / B(2); 
    
    [rocx,rocy,roct,auc] = perfcurve(Y,scores,1);
    
    hs(e)=subplot(2,3,e); hold on
    %bins=linspace(0,1.25,100);
    bins=linspace(-13,0.25,100);
    [N1,X1]=hist(absIn,bins);
    [N2,X2]=hist(absOut,bins);
    hold on
    relfreq1=N1/sum(N1);
    relfreq2=N2/sum(N2);
    %thresh10=X1(find(relfreq1./relfreq2>10,1));
    hp(e,1)=area(X2,relfreq2,'FaceAlpha',0.5);
    hp(e,2)=area(X1,relfreq1,'FaceAlpha',0.5);
    hthr=plot(thresh*[1 1],[0 0.125],'--k');
    hthr10=plot(thresh95*[1 1],[0 0.1],':k');
    axis square
    
    htxt(e)=text(thresh+0.05,0.125,[num2str(exp(thresh),'%0.2f') ' J/cm^3'] );
    htxt10(e)=text(thresh95+0.05,0.1,[num2str(exp(thresh95),'%0.2f') ' J/cm^3'] );

    htit(e)=title(sprintf('Echo %d',e),'FontWeight','normal','FontSize',16);
    
    hs(3+e)=subplot(2,3,3+e); hold on
    plot(rocx,rocy);
    plot([0 1],[0 1],'--k');
    htxt(e)=text(0.2,0.8,sprintf(['AUC=%0.2f'],auc));
    axis square
    
    if e==1
        hlg=legend([hp(e,1) hp(e,2) hthr hthr10],'n.s.','sig.','p=0.5','p=0.95');
    end
    
    
    
end

%% esthetics and print
set(htxt(1:3),'FontSize',12);

set(hs(1:3),'YLim',[0 0.2]);
set(hlg,'box','off');
set(get(hs(2),'XLabel'),'String','Log Absorption (J/cm^3)','FontSize',fs);
set(get(hs(1),'YLabel'),'String','Rel. Freq.','FontSize',fs);
set(get(hs(5),'XLabel'),'String','False Positive Rate','FontSize',fs);
set(get(hs(4),'YLabel'),'String','True Positive Rate','FontSize',fs);


%set(hs(1:3),'XTick',[0:0.25:1.25]);
set(hs(1:3),'YTick',[0:0.05:0.2]);
set(hs(2),'YTickLabel',{'','','','',''});
set(hs(3),'YTickLabel',{'','','','',''});


set(hs(4:6),'XTick',[0:0.2:1]);
set(hs(4:6),'YTick',[0:0.2:1]);
set(hs(5),'YTickLabel',{'','','','',''});
set(hs(6),'YTickLabel',{'','','','',''});

delUp=0;
delLeft=0.04;
for s=2:3:6
    subPos=get(hs(s),'Position');
    set(hs(s),'Position',[subPos(1)-delLeft subPos(2) subPos(3) subPos(4)]);
end
for s=3:3:6
    subPos=get(hs(s),'Position');
    set(hs(s),'Position',[subPos(1)-2*delLeft subPos(2) subPos(3) subPos(4)]);
end
for s=4:6
    subPos=get(hs(s),'Position');
    set(hs(s),'Position',[subPos(1) subPos(2)+delUp subPos(3) subPos(4)]);
end
%%

%%
lgPos=get(hlg,'Position');
set(hlg,'Position',[0.13 0.795 lgPos(3) lgPos(4)]);

%%
for e=1:3
    titPos=get(htit(e),'Position');
    set(htit(e),'Position',[titPos(1) titPos(2)+0.01 titPos(3)]);
    set(htit(e),'FontWeight','normal');
end

%%
sublabel([hs(1:3:6)],-10,-30,'FontSize',16,'FontWeight','bold');
print -dpng ../figures/thresholdAnalysis
crop('../figures/thresholdAnalysis.png',0);







