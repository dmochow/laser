clear all; close all; clc

%% coordinates of average laser origin
contourLevels=[0.01 0.1 1 10];
stLaserOrigin = round([56.7901 4.8800 87.3221]);
stUlf=[64.2489,-2.1050,97.0123];
stUrf=[44.7676,2.9479,96.1140];
stLaserOrigin(1)=161-stLaserOrigin(1);
stUlf(1)=161-stUlf(1);
stUrf(1)=161-stUrf(1);
colors=[0 0.4470 0.7410;0.8500    0.3250    0.0980];
%% colormap
% cmap = cbrewer('div','PuOr',100);
% cm=flipud(cmap);
% cm2=cat(1,[1 1 1],cm);

cm=colormap('jet');

%% anatomy
pathToData='../data/SAVG/NII';
anatFilename=fullfile(pathToData,'TT_N27+tlrc');
brainMaskFilename=fullfile(pathToData,'brain_mask+tlrc.BRIK');
absorptionFilename=fullfile(pathToData,'Amu+tlrc.BRIK');
[~, anat, info, ~] = BrikLoad (anatFilename);
anat=anat(end:-1:1,:,:);
[~, brainMask, ~, ~] = BrikLoad (brainMaskFilename);
brainMask=brainMask(end:-1:1,:,:);
[~, absorption, ~, ~] = BrikLoad (absorptionFilename);
absorption=absorption(end:-1:1,:,:);
[~,labels,~,~]=BrikLoad('../data/SAVG/NII/Classes+tlrc.BRIK');
labels=labels(end:-1:1,:,:);


slice2show=88;
laserOrigin=[104 5];
mrdims=size(anat);
anat(~brainMask)=NaN;
absorption(~brainMask)=NaN;
logAbsorption=log10(absorption);
%logAbsorption(logAbsorption==-Inf)=NaN;

%%
hf=figure;
hs(1)=subplot(221);
hold on
xmri=-(mrdims(1)-1)/2:(mrdims(1)-1)/2;
ymri=-(mrdims(2)-1)/2:(mrdims(2)-1)/2;
[XMRI,YMRI]=ndgrid(1:mrdims(1),1:mrdims(2));
hpc=pcolor(XMRI,YMRI,anat(:,:,slice2show));
set(hpc,'EdgeColor','none');
%set(hpc,'FaceAlpha',0.5);
axis ij;
axis equal;
%axis off
colormap bone;
freezeColors;
%logAbsorption(logAbsorption(:,:,slice2show)<0.001)=NaN;
hpc=pcolor(XMRI,YMRI,logAbsorption(:,:,slice2show));
set(hpc,'EdgeColor','none');
set(hpc,'FaceAlpha',0.65);
colormap(cm);
caxis([-3 0]);
hcb=colorbar;
set(hcb,'Ticks',[-4:1:0]);
set(hcb,'TickLabels',{'0.0001','0.001','0.01','0.1','1'});
hyl=ylabel(hcb,'J/cm^3');
set(hyl,'Rotation',0);
ylPos=get(hyl,'Position');
set(hyl,'Position',[0.6 0.65 0]);
set(gca,'XTickLabelMode','auto');
set(gca,'YTickLabelMode','auto');
xlim([15 150]);
ylim([0 180]);
%hxl=title(gca,'z=88');set(hxl,'Visible','on','FontWeight','normal')
axis off
%axPos=get(gca,'Position');
%set(hs(1),'Position',[axPos(1) axPos(2)+0.05 axPos(3) axPos(4)]);

lgPos=get(hyl,'Position');
set(hyl,'Position',[lgPos(1) lgPos(2)-0.3 lgPos(3)]);

hpl=plot([stUrf(1) stUlf(1)],[stUrf(2) stUlf(2)],'LineWidth',3,'Color',[1 0 0]);
htxt=text(stUrf(1)-40,stUrf(2)-12,'Illumination');

%%
RFP=[95 22];
[~,xcoor]=min(abs(XMRI(:,1)-RFP(1)));
[~,ycoor]=min(abs(YMRI(1,:)-RFP(2)));
A_rfp=absorption(xcoor,ycoor,slice2show);

RDLPFC=[123 36];
[~,xcoor]=min(abs(XMRI(:,1)-RDLPFC(1)));
[~,ycoor]=min(abs(YMRI(1,:)-RDLPFC(2)));
A_rdlpfc=absorption(xcoor,ycoor,slice2show);

LDLPFC=[46 33];
[~,xcoor]=min(abs(XMRI(:,1)-LDLPFC(1)));
[~,ycoor]=min(abs(YMRI(1,:)-LDLPFC(2)));
A_ldlpfc=absorption(xcoor,ycoor,slice2show);

LFP=[79 22];
[~,xcoor]=min(abs(XMRI(:,1)-LFP(1)));
[~,ycoor]=min(abs(YMRI(1,:)-LFP(2)));
A_lfp=absorption(xcoor,ycoor,slice2show);

MPFC=[84 30];
[~,xcoor]=min(abs(XMRI(:,1)-MPFC(1)));
[~,ycoor]=min(abs(YMRI(1,:)-MPFC(2)));
A_mpfc=absorption(xcoor,ycoor,slice2show);

hs(2)=subplot(222);
bar(1:4,[A_rfp A_rdlpfc A_mpfc A_lfp]);
set(gca,'XTickLabel',{'R-Frontal Pole','R-DLPFC','mPFC','L-Frontal Pole','L-DLPFC'});
set(gca,'XTickLabelRotation',45)
ylabel('Absorption (J/cm^3)');
box off

sublabel(hs,-20,0,'FontWeight','bold');
print -dpng ../figures/figureDan
crop('../figures/figureDan.png',0);