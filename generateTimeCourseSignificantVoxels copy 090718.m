clear all; close all; clc
e=2;
nTRsToSkip=5; % while gradients stabilize
TR=2.8;
time=(0:644)*TR;
pathToData='../data/SAVG/NII';
anatFilename=fullfile(pathToData,'TT_N27+tlrc');
muBoldFilename=fullfile(pathToData,['s8muwBoldEcho' num2str(e) '-WhiteMatter1-KWHM3-CSF0-KCSF1-Center1+tlrc.BRIK']);
brainMaskFilename=fullfile(pathToData,'brain_mask+tlrc.BRIK');
isSigBefDurFilename=fullfile(pathToData,['thresh_isSigBefDur_FWHM8_ECHO'...
    num2str(e) '-WhiteMatter1-KWM3-CSF0-KCSF1-Center1+tlrc.BRIK']);
isSigBefAftFilename=fullfile(pathToData,['thresh_isSigBefAft_FWHM8_ECHO'...
    num2str(e) '-WhiteMatter1-KWM3-CSF0-KCSF1-Center1+tlrc.BRIK']);

[~,muBold,~,~]=BrikLoad(muBoldFilename);

% throw out first 5 TRs
muBold(:,:,:,1:nTRsToSkip)=0;

[~,isSigBefDur,~,~]=BrikLoad(isSigBefDurFilename);
[~,isSigBefAft,~,~]=BrikLoad(isSigBefAftFilename);

sigTs1=vol2ts(muBold,isSigBefDur>0);
sigTs2=vol2ts(muBold,isSigBefAft>0);

muTs1=mean(sigTs1,2);
muTs2=mean(sigTs2,2);

%f = fit(time.',mean(sigTs1,2),'smoothingspline','SmoothingParam',0.01);

winsize=20;
B1 = smoothdata(muTs1,1,'rloess',winsize);
B2 = smoothdata(muTs2,1,'rloess',winsize);

%%
cmap = cbrewer('div','PuOr',100);
cm=colormap(flipud(cmap));

%%
hf=figure;  
%subplot(211);
%subplot(221); plot(f,time.',mean(sigTs1,2)); %plot(mean(sigTs1,2));
hold on 
hraw=plot(time,muTs1,'.','LineStyle','none','MarkerEdgeColor',[0.7 0.7 0.7]); %plot(time,muTs1);
%plot(time,muTs1,'.','LineStyle','none','MarkerEdgeColor',cm(end,:)); %plot(time,muTs1);
hsm=plot(time,B1,'LineWidth',2,'Color',[1 0 0]);
ylim([-0.1 0.1]);
yl=ylim;
set(gca,'YTick',[-0.1 0 0.1]);
harea=area([600 1200],[yl(2) yl(2)],'BaseValue',yl(1),'FaceColor',[0.7 0.7 0.7],'EdgeColor','none','FaceAlpha',0.25);
xlabel('Time (mm:ss)','FontSize',20);
ylabel('BOLD (a.u.)','FontSize',20);
set(gca,'XTick',[0:300:1800]);
set(gca,'XTickLabel',{'0:00','5:00','10:00','15:00','20:00','25:00','30:00'});
xlim([0 1800]);
hlg=legend([hsm harea],'smoothed','laser on');
set(hlg,'box','off');
set(hlg,'FontSize',20);
print -dpng ../figures/ts1.png
crop('../figures/ts1.png',0);

%%
hf=figure;  
hold on 
hraw=plot(time,muTs2,'.','LineStyle','none','MarkerEdgeColor',[0.7 0.7 0.7]); %plot(time,muTs1);
hsm=plot(time,B2,'LineWidth',2,'Color',[1 0 0]);
ylim([-0.1 0.1]);
yl=ylim;
set(gca,'YTick',[-0.1 0 0.1]);
harea=area([600 1200],[yl(2) yl(2)],'BaseValue',yl(1),'FaceColor',[0.7 0.7 0.7],'EdgeColor','none','FaceAlpha',0.25);
xlabel('Time (mm:ss)','FontSize',20);
ylabel('BOLD (a.u.)','FontSize',20);
set(gca,'XTick',[0:300:1800]);
set(gca,'XTickLabel',{'0:00','5:00','10:00','15:00','20:00','25:00','30:00'});
xlim([0 1800]);
hlg=legend([hsm harea],'smoothed','laser on');
set(hlg,'box','off');
set(hlg,'FontSize',20);
print -dpng ../figures/ts2.png
crop('../figures/ts2.png',0);