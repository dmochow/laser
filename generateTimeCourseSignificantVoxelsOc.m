% 10.19.19 look at time course of optimally-combined bold

clear all; close all; clc
PRINT_FIGURE=0;
nTRs=645;
LASERONSETTR=215;
LASEROFFSETTR=430;
testInds=setdiff(1:nTRs,[1:LASERONSETTR  LASEROFFSETTR:nTRs]);
nPerms=1000;
nTRsToSkip=5; % while gradients stabilize
TR=2.8;
time=(0:nTRs-1)*TR;
pathToData='../data/SAVG/NII';
anatFilename=fullfile(pathToData,'TT_N27+tlrc');
muBoldFilenames={'s8muocBold-Filter1-Derivative1-WhiteMatter3-STANDARDIZE1-TSHIFT1-CSHIFT2+tlrc.BRIK'};
filenames={'thresh_oc_isSigBefDur_NLAGS5_CONSTANT0_s8muocBold-Filter1-Derivative1-WhiteMatter3-STANDARDIZE1-TSHIFT1-CSHIFT2+tlrc.BRIK'};
arFilenames={'oc_arCoeffsResidsEcho-NLAGS5_CONSTANT0.mat'};
coeffIndsDim3=[1;1;1;2]; % 1=pre vs stim; 2=pre vs post
brainMaskFilename='resampled_brain_mask+tlrc.BRIK';
lineColor=[0 0 0];
winsize=20; % time series smoothing parameter
legDelRight=0.075;
legDelUp=0.025;
nFilenames=numel(filenames);
laserOrigin = [40 1 34]+[1 1 1]; % LR, AP, IS
nfft=32; % for AR spectra
fs=1/TR; % for AR spectra
colors=[0 0.4470 0.7410;0.8500    0.3250    0.0980];

[~,brainMask,~,~]=BrikLoad(fullfile(pathToData,brainMaskFilename));
brainMask=logical(brainMask);



hf=figure;
for e=1:nFilenames
    [~,muBold,~,~]=BrikLoad(fullfile(pathToData,muBoldFilenames{e}));
    muBold(:,:,:,1:nTRsToSkip)=0;
    
    thisBefDurFilename=['thresh_oc_isSigBefDur_NLAGS5_CONSTANT0_s8muocBold-Filter1-Derivative1-WhiteMatter3-STANDARDIZE1-TSHIFT1-CSHIFT2+tlrc.BRIK'];
    thisBefAftFilename=['thresh_oc_isSigBefAft_NLAGS5_CONSTANT0_s8muocBold-Filter1-Derivative1-WhiteMatter3-STANDARDIZE1-TSHIFT1-CSHIFT2+tlrc.BRIK'];
    
    [~,isSigBefDur,~,~]=BrikLoad(fullfile(pathToData,thisBefDurFilename));
    sigBefDurTs=vol2ts(muBold,isSigBefDur>0);
    muSigBefDurTs=mean(sigBefDurTs,2);
    
    [~,isSigBefAft,~,~]=BrikLoad(fullfile(pathToData,thisBefAftFilename));
    sigBefAftTs=vol2ts(muBold,isSigBefAft>0);
    muSigBefAftTs=mean(sigBefAftTs,2);
    
    
    % AR coefficients
    load(fullfile(pathToData,arFilenames{e}),'coeffsPreStim','coeffsPrePost');
    
    freqs=(0:nfft-1)/nfft*fs;
    isSig1DbefDur=vol2ts(isSigBefDur,brainMask);
%     isSig1DbefAft=vol2ts(isSigBefAft,brainMask);
    
    % before vs during
    coeffsPreStim=-coeffsPreStim; % 
    coeffsPreStim=cat(1,ones(1,size(coeffsPreStim,2),size(coeffsPreStim,3)),coeffsPreStim);
    A=fft(coeffsPreStim,nfft,1);
    C=1./(A+eps); %invert to get the AR spectrum
    absC=abs(C);
    tmpSbefDur=squeeze(mean(absC(:,isSig1DbefDur==1,:),2));
    tmpNbefDur=squeeze(mean(absC(:,isSig1DbefDur==0,:),2));
    
    % before vs after
    coeffsPrePost=-coeffsPrePost; % 
    coeffsPrePost=cat(1,ones(1,size(coeffsPrePost,2),size(coeffsPrePost,3)),coeffsPrePost);
    A=fft(coeffsPrePost,nfft,1);
    C=1./(A+eps); %invert to get the AR spectrum
    absC=abs(C);
%     tmpSbefAft=squeeze(mean(absC(:,isSig1DbefAft==1,:),2));
%     tmpNbefAft=squeeze(mean(absC(:,isSig1DbefAft==0,:),2));
    
    % 09/21/18: substitute which voxels we shoe
    tmpSbefAft=squeeze(mean(absC(:,isSig1DbefDur==1,:),2));
    tmpNbefAft=squeeze(mean(absC(:,isSig1DbefDur==0,:),2));
    
    
    % 10/17/18: compute percent signal change
    %%
    for v=1:size(sigBefDurTs,2)
        fBold{e}(:,v)=filter(1,coeffsPreStim(:,v,2),sigBefDurTs(:,v));
    end
    mufBold(:,e)=mean(fBold{e},2);
    
%     for v=1:size(sigBefDurTs,2)
%         fBold{e,2}(:,v)=filter(1,coeffsPrePost(:,v,2),sigBefAftTs(:,v));
%     end
    %%
    
    
    % row 1: time series
    hs(e)=subplot(4,3,e);hold on
    % show AR filtered
%     hp1=plot(time(1:LASERONSETTR),mufBold(1:LASERONSETTR,e),'Color',colors(1,:));
%     hp2=plot(time(LASERONSETTR+1:LASEROFFSETTR),mufBold(LASERONSETTR+1:LASEROFFSETTR,e),'Color',colors(2,:));
%     hp3=plot(time(LASEROFFSETTR+1:end),mufBold(LASEROFFSETTR+1:end,e),'Color',[0.7 0.7 0.7]);

    % show the BOLD
    hp1=plot(time(1:LASERONSETTR),muSigBefDurTs(1:LASERONSETTR),'Color',colors(1,:));
    hp2=plot(time(LASERONSETTR+1:LASEROFFSETTR),muSigBefDurTs(LASERONSETTR+1:LASEROFFSETTR),'Color',colors(2,:));
    hp3=plot(time(LASEROFFSETTR+1:end),muSigBefDurTs(LASEROFFSETTR+1:end),'Color',[0.7 0.7 0.7]);
    if e==1
        hlg1=legend('Pre','tPBM','Post');
        set(hlg1,'box','off');
    end
    %axis square
    
    % row 2: AR spectra of sig voxels
    hs(e+3)=subplot(4,3,e+3);
    hstem=stem(freqs,[tmpSbefDur(:,1) tmpSbefDur(:,2) tmpSbefAft(:,2)],'filled');
    set(hstem(3),'MarkerEdgeColor',[0.7 0.7 0.7]);
    set(hstem(3),'MarkerFaceColor',[0.7 0.7 0.7]);
    %axis square
    xlim([0 fs/2])
    if e==1
        hlg2=legend('Pre','tPBM','Post');
        set(hlg2,'box','off');
    end
    
    % 10.16.18 percent signal change
    posFreqMask=1:nfft/2+1;
    freqWeighting=tmpSbefDur(1:nfft/2+1,1)/sum(tmpSbefDur(1:nfft/2+1,1));
    pscf=( tmpSbefDur(1:nfft/2+1,2)-tmpSbefDur(1:nfft/2+1,1) )./tmpSbefDur(1:nfft/2+1,1);
    psc(e)=sum(freqWeighting .* pscf);
    pscf2=( tmpSbefAft(1:nfft/2+1,2)-tmpSbefAft(1:nfft/2+1,1) )./tmpSbefAft(1:nfft/2+1,1);
    psc2(e)=sum(freqWeighting .* pscf2);
    
    % row 3: AR spectra of n.s. voxels
    hs(e+6)=subplot(4,3,e+6);
    hstem=stem(freqs,[tmpNbefDur(:,1) tmpNbefDur(:,2) tmpNbefAft(:,2)],'filled');
    set(hstem(3),'MarkerEdgeColor',[0.7 0.7 0.7]);
    set(hstem(3),'MarkerFaceColor',[0.7 0.7 0.7]);
    xlim([0 fs/2])
    
    hs(e+9)=subplot(4,3,e+9);
    hbar(e)=bar([1 2],[psc(e) psc2(e)]*100);
    
end

% % axis labels
% set(get(hs(1),'ylabel'),'String','BOLD (a.u.)');
% set(get(hs(4),'ylabel'),'String','AR Spectrum');
% set(get(hs(7),'ylabel'),'String','AR Spectrum');
% set(get(hs(10),'ylabel'),'String','Percent Change');
% 
% set(get(hs(1),'xlabel'),'String','Time (mm:ss)');
% set(get(hs(2),'xlabel'),'String','Time (mm:ss)');
% set(get(hs(3),'xlabel'),'String','Time (mm:ss)');
% set(get(hs(7),'xlabel'),'String','Frequency (Hz)');
% set(get(hs(8),'xlabel'),'String','Frequency (Hz)');
% set(get(hs(9),'xlabel'),'String','Frequency (Hz)');
% 
% % axis limits
% set(hs(1:3),'XLim',[0 1800]);
% set(hs(1:3),'YLim',[-0.2 0.2]);
% set(hs(4:9),'YLim',[0 6]);
% set(hs(10:12),'YLim',[0 40]);
% % set(hs(4:9),'XLim',[0 0.1]); % uncomment if you want to show only
% % passband
% %set(hs(4:9),'YLim',[0 6]);
% 
% % x-axis ticks
% set(hs(1:3),'XTick',[600:600:1800]);
% set(hs(1:3),'XTickLabel',{'10:00','20:00','30:00'});
% set(hs(10:12),'XTick',[1 2]);
% set(hs(10:12),'XTickLabel',{'tPBM','Post'});
% 
% 
% % y-axis ticks
% set(hs([2 3 5 6 8 9]),'YTickLabel',{'','','','',''});
% set(hs([11:12]),'YTickLabel',{'','',''});
% 
% % titles
% set(get(hs(1),'Title'),'String','Echo 1');
% set(get(hs(2),'Title'),'String','Echo 2');
% set(get(hs(3),'Title'),'String','Echo 3');
% set(get(hs(1),'Title'),'FontWeight','normal');
% set(get(hs(2),'Title'),'FontWeight','normal');
% set(get(hs(3),'Title'),'FontWeight','normal');
% 
% % boxes
% set(hs(:),'Box','off');
% 
% % subplot positions
% xMargin=0.1;
% xStep=(1-xMargin)/3;
% subposWidth=0.275;
% 
% for s=1:3:9
%     subPos=get(hs(s),'Position');
%     set(hs(s),'Position',[xMargin subPos(2) subposWidth subPos(4)]);
% end
% for s=2:3:9
%     subPos=get(hs(s),'Position');
%     set(hs(s),'Position',[xMargin+xStep subPos(2) subposWidth subPos(4)]);
% end
% for s=3:3:9
%     subPos=get(hs(s),'Position');
%     set(hs(s),'Position',[xMargin+2*xStep subPos(2) subposWidth subPos(4)]);
% end
% 
% delUp=0.025;
% for s=1:3
%     subPos=get(hs(s),'Position');
%     set(hs(s),'Position',[subPos(1) subPos(2)+delUp subPos(3) subPos(4)]);
% end
% 
% delDown=0.025;
% for s=4:6
%     subPos=get(hs(s),'Position');
%     set(hs(s),'Position',[subPos(1) subPos(2)-delDown subPos(3) subPos(4)]);
% end
% 
% delDown=0.05;
% for s=10:12
%     subPos=get(hs(s),'Position');
%     set(hs(s),'Position',[subPos(1) subPos(2)-delDown subPos(3) subPos(4)]);
% end
% 
% % move legend up
% lgDelUp=0.03;
% lgPos=get(hlg1,'Position');
% set(hlg1,'Position',[lgPos(1)-0.14 lgPos(2)+0.05 lgPos(3) lgPos(4)]);
% 
% lgPos=get(hlg2,'Position');
% set(hlg2,'Position',[lgPos(1) lgPos(2)+lgDelUp lgPos(3) lgPos(4)]);
% 
% % remove edges from bar
% set(hbar(:),'EdgeColor','none');
% %
% %%
% if PRINT_FIGURE
% sublabel([hs(1) hs(4) hs(7) hs(10)],-10,-45,'FontSize',16,'FontWeight','Bold');
% print -dpng ../figures/ocTimeSeriesFigure.png
% end
% crop('../figures/ocTimeSeriesFigure.png',0);