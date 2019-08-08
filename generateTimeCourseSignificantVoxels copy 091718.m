clear all; close all; clc

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
muBoldFilenames={'s8muwBoldEcho1-WhiteMatter1-KWHM3-CSF0-KCSF1-Center1+tlrc.BRIK',...
    's8muwBoldEcho2-WhiteMatter1-KWHM3-CSF0-KCSF1-Center1+tlrc.BRIK',...
    's8muwBoldEcho3-WhiteMatter1-KWHM3-CSF0-KCSF1-Center1+tlrc.BRIK',...
    's8muwBoldEcho1-WhiteMatter1-KWHM3-CSF0-KCSF1-Center1+tlrc.BRIK'};
filenames={'thresh_isSigBefDur_FWHM8_ECHO1-WhiteMatter1-KWM3-CSF0-KCSF1-Center1+tlrc.BRIK',...
    'thresh_isSigBefDur_FWHM8_ECHO2-WhiteMatter1-KWM3-CSF0-KCSF1-Center1+tlrc.BRIK',...
    'thresh_isSigBefDur_FWHM8_ECHO3-WhiteMatter1-KWM3-CSF0-KCSF1-Center1+tlrc.BRIK',...
    'thresh_isSigBefAft_FWHM8_ECHO1-WhiteMatter1-KWM3-CSF0-KCSF1-Center1+tlrc.BRIK'};
arFilenames={'arCoeffsResidsEcho1.mat','arCoeffsResidsEcho2.mat','arCoeffsResidsEcho3.mat','arCoeffsResidsEcho1.mat'};
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
for f=1:nFilenames
    [~,muBold,~,~]=BrikLoad(fullfile(pathToData,muBoldFilenames{f}));
    muBold(:,:,:,1:nTRsToSkip)=0;
    
    [~,isSig,~,~]=BrikLoad(fullfile(pathToData,filenames{f}));
    sigTs=vol2ts(muBold,isSig>0);
    muSigTs=mean(sigTs,2);
    smuSigTs = smoothdata(muSigTs,1,'rloess',winsize);

    % AR coefficients
    load(fullfile(pathToData,arFilenames{f}),'coeffsPreStim','coeffsPrePost','rsdPreStim','rsdPrePost');
    
    if coeffIndsDim3(f)==1
        coeffs=coeffsPreStim;
    else
        coeffs=coeffsPrePost;
    end
    % 

    freqs=(0:nfft-1)/nfft*fs;
    isSig1D=vol2ts(isSig,brainMask);

    coeffs(1,:,:)=1; % hack because idk what I'm doing yet
    coeffs(2:end,:,:)=-coeffs(2:end,:,:); % 
    A=fft(coeffs,nfft,1);
    C=1./(A+eps); %invert to get the AR spectrum
    absC=abs(C);
    
    tmpS=squeeze(mean(absC(:,isSig1D==1,:),2));
    tmpN=squeeze(mean(absC(:,isSig1D==0,:),2));
    
    rsdS=squeeze(mean(rsdPreStim(:,isSig1D==1,:),2));
    rsdN=squeeze(mean(rsdPreStim(:,isSig1D==0,:),2));
    
    hs(f+4)=subplot(2,4,f+4);
    hstem=stem(freqs,[tmpS(:,1) tmpS(:,2)],'filled');
    axis square
    xlim([0 fs/2])
    if coeffIndsDim3(f)==1
        hlg(f+4)=legend('Pre','Stim');
    else
        hlg(f+4)=legend('Pre','Post');
    end
    set(hlg(f+4),'box','off');
    
    hs(f)=subplot(2,4,f);hold on
    if coeffIndsDim3(f)==1
        hp1=plot(time(1:LASERONSETTR),muSigTs(1:LASERONSETTR),'Color',colors(1,:));
        hp2=plot(time(LASERONSETTR+1:LASEROFFSETTR),muSigTs(LASERONSETTR+1:LASEROFFSETTR),'Color',colors(2,:));
        hp3=plot(time(LASEROFFSETTR+1:end),muSigTs(LASEROFFSETTR+1:end),'Color',[0.7 0.7 0.7]);
        hlg(f)=legend('Pre','Stim');
        htit(f)=title(['Echo ' num2str(f)],'FontWeight','normal');
        axis square
    else
        hp1=plot(time(1:LASERONSETTR),muSigTs(1:LASERONSETTR),'Color',colors(1,:));
        hp3=plot(time(LASEROFFSETTR+1:end),muSigTs(LASEROFFSETTR+1:end),'Color',colors(2,:));
        hp2=plot(time(LASERONSETTR+1:LASEROFFSETTR),muSigTs(LASERONSETTR+1:LASEROFFSETTR),'Color',[0.7 0.7 0.7]);
        hlg(f)=legend('Pre','Post');
        htit(f)=title(['Echo 1' ],'FontWeight','normal');
        axis square
    end
    set(hlg(f),'box','off');
    
    
    %plot([rsdS(:,1)]);
    %plot([rsdS(:,1) rsdN(:,1)]);
    
    % show as supplementary figure
%     subplot(2,4,f+4);
%     plot(freqs,[tmpN(:,1) tmpN(:,2)]);
    
    
%     subplot(4,2,f*2);
%     stem(freqs,[tmpS(:,2) tmpN(:,2)]);
    
%     stem(freqs,[tmpS(:,1) tmpN(:,1)]);
%     legend('Sig Voxels','N.S. Voxels');
%     subplot(4,2,f*2);
%     stem(freqs,[tmpS(:,2) tmpN(:,2)]);
    
    
%     % permutation testing
%     mockSmuSigTs=zeros(size(muSigTs,1),nPerms);
%     for p=1:nPerms
%         p
%         tmp=surrogateResponseGenerator(muSigTs);
%         mockSmuSigTs(:,p)= smoothdata(tmp,1,'rloess',winsize);
%     end
%     
%     [mockSems,mockMus,mockStds]=nansem(mockSmuSigTs,2);
%         
%     hs(f)=subplot(3,2,f);
%     hold on
%     
%     hplot(f)=plot(time,smuSigTs,'LineWidth',2,'Color',lineColor);
%     hmock(f)=shadedErrorBar(time,mockMus,1.96*mockStds,'k',1);
%     set(hmock(f).mainLine,'LineStyle','none');
%     set(hmock(f).patch,'FaceColor',[0,0,0.5]);
%     ylim([-0.1 0.1]);
%     yl=ylim;
%     set(gca,'YTick',[-0.1 0 0.1]);
%     harea=area([600 1200],[yl(2) yl(2)],'BaseValue',yl(1),'FaceColor',[0.7 0.7 0.7],'EdgeColor','none','FaceAlpha',0.25);
% 
%     
%     set(gca,'XTick',[0:600:1800]);
%     set(gca,'XTickLabel',{'0:00','10:00','20:00','30:00'});
%     xlim([0 1800]);
%     %axis square
%     if f==1
%     hlg=legend([hplot(1) hmock(1).patch],'Sig. Voxels','95% C.I.');
%     set(hlg,'box','off');
%     end
%     %set(hlg,'FontSize',14);
     
end

set(get(hs(1),'ylabel'),'String','BOLD (a.u.)');
set(get(hs(5),'ylabel'),'String','Magnitude Spectrum');
% set(get(hs(3),'ylabel'),'String','BOLD (a.u.)');
% set(get(hs(3),'ylabel'),'FontSize',20);
% 


set(get(hs(5),'xlabel'),'String','Time (s)');
set(get(hs(6),'xlabel'),'String','Time (s)');
set(get(hs(7),'xlabel'),'String','Time (s)');
set(get(hs(8),'xlabel'),'String','Time (s)');

set(get(hs(5),'xlabel'),'String','Frequency (Hz)');
set(get(hs(6),'xlabel'),'String','Frequency (Hz)');
set(get(hs(7),'xlabel'),'String','Frequency (Hz)');
set(get(hs(8),'xlabel'),'String','Frequency (Hz)');
% 
% 
% set(get(hs(1),'Title'),'String','Echo 1');
% set(get(hs(2),'Title'),'String','Echo 2');
% set(get(hs(3),'Title'),'String','Echo 3');
% set(get(hs(4),'Title'),'String','Echo 3, Pre vs Post');
% set(get(hs(1),'Title'),'FontWeight','normal');
% set(get(hs(2),'Title'),'FontWeight','normal');
% set(get(hs(3),'Title'),'FontWeight','normal');
% set(get(hs(4),'Title'),'FontWeight','normal');
% 
% heightScale=1.15;
% for s=1:4
%     subPos=get(hs(s),'Position');
%     set(hs(s),'Position',[subPos(1) subPos(2) subPos(3) subPos(4)*heightScale]);
% end
% 
delUp=0.2;
for s=5:8
    subPos=get(hs(s),'Position');
    set(hs(s),'Position',[subPos(1) subPos(2)+delUp subPos(3) subPos(4)]);
end
% 
% % move legend up
% lgDelUp=0.02;
%    lgPos=get(hlg,'Position');
%     set(hlg,'Position',[lgPos(1) lgPos(2)+lgDelUp lgPos(3) lgPos(4)]);
% 
%%
sublabel([hs(1) hs(5)],-10,-30,'FontSize',16,'FontWeight','Bold');
print -dpng ../figures/timeSeriesFigure.png
crop('../figures/timeSeriesFigure.png',0);