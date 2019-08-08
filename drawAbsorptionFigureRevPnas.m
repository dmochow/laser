clear all; close all; clc

%% coordinates of average laser origin
contourLevels=[0.01 0.1 1 10];

%% colormap
cmap = cbrewer('div','PuOr',100);
cm=flipud(cmap);
cm2=cat(1,[1 1 1],cm);

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

%% load absorption error bars
load('../data/precomputed/absorptionErrorBars','csfTotals','gmTotals','wmTotals','gmMaxima','wmMaxima','csfMaxima');
[semGm,muGm]=nansem(gmTotals);
[semWm,muWm]=nansem(wmTotals);
[semCsf,muCsf]=nansem(csfTotals);


slice2show=88;
laserOrigin=[104 5];
mrdims=size(anat);
anat(~brainMask)=NaN;
absorption(~brainMask)=NaN;
logAbsorption=log10(absorption);
logAbsorption(logAbsorption==-Inf)=NaN;

%%
hf=figure;
ha = tight_subplot(2,2,[.05 .05],[.1 .1],[.1 .1])
%for ii = 1:6; axes(ha(ii)); plot(randn(10,ii)); end
%set(ha(1:4),'XTickLabel',''); set(ha,'YTickLabel','')
%axis off
axes(ha(1));
axis off
%set(ha(1),'YAxisLocation','right');
%set(ha(1),'YAxisLocation','origin');
%set(ha(1),'XAxisLocation','top');
hold on
% add sidebar showing layers
SCALP=9.45;
SKULL=8.8;
CSF=8.2;
BRAIN=8;
% har1=area([1 2],[SCALP SCALP],'BaseValue',SKULL,'FaceAlpha',0.5,'LineStyle','none');
% har2=area([1 2],[SKULL SKULL],'BaseValue',CSF,'FaceAlpha',0.5,'LineStyle','none');
% har3=area([1 2],[CSF CSF],'BaseValue',BRAIN,'FaceAlpha',0.5,'LineStyle','none');
% har4=area([1 2],[BRAIN BRAIN],'BaseValue',0,'FaceAlpha',0.5,'LineStyle','none');
har1=area([-1 1],[SCALP SCALP],'BaseValue',SKULL,'FaceAlpha',0.5,'LineStyle','none');
har2=area([-1 1],[SKULL SKULL],'BaseValue',CSF,'FaceAlpha',0.5,'LineStyle','none');
har3=area([-1 1],[CSF CSF],'BaseValue',BRAIN,'FaceAlpha',0.5,'LineStyle','none');
har4=area([-1 1],[BRAIN BRAIN],'BaseValue',0,'FaceAlpha',0.5,'LineStyle','none');

xcoords=[0.16:0.01:0.18];
HT=0.83;
for p=1:3
    if p==2
        hann(p)=annotation('textarrow',[xcoords(p) xcoords(p)],[0.06+HT HT],'String','Photons');
        set(hann(p),'Color','r');
    else
        hann(p)=annotation('textarrow',[xcoords(p) xcoords(p)],[0.06+HT HT],'String','');
        set(hann(p),'Color','r');
    end
end

hann=annotation('textarrow',[0.17 0.34],[HT HT],'String','');
set(hann,'Color','k');
hann=annotation('textarrow',[0.17 0.17],[HT 0.627],'String','');
set(hann,'Color','k');

htext=text(0.06,1,'Depth','FontSize',12);
htext=text(0.5,10.45,'Radius','FontSize',12);

% hann=annotation('textarrow',[0.31 0.31],[0.94 0.88],'String','');
% set(hann,'Color','r');
% hann=annotation('textarrow',[0.27 0.27],[0.94 0.88],'String','');
% set(hann,'Color','r');

%annotation('textarrow',[0.30 0.30],[0.965 0.93],'String','');
%annotation('textarrow',[0.28 0.28],[0.965 0.93],'String','');
set(gca,'XLimMode','auto');
set(gca,'YLimMode','auto');
%xlim([1 2]);
ylim([0.01 SCALP]);

set(gca,'XTick',[0.25]);
set(gca,'XTickLabel',{'20 mm'});

% set(gca,'YTick',[BRAIN,SCALP]);
% set(gca,'YTickLabel',{'80 mm','94.5'});

set(gca,'YTick',[0 7.45]);
set(gca,'YTickLabel',{'94.5','20 mm'});

axPos=get(ha(1),'Position');
set(ha(1),'Position',[axPos(1) axPos(2)+0.1 axPos(3) axPos(4)/1.5]);
x = [0 eps eps 0];
y = [0 0 eps eps];
hlg0=legend([har1 har2 har3 har4],'Scalp','Skull','CSF','Brain');
set(hlg0,'orientation','horizontal');
set(hlg0,'box','off');
lgPos=get(hlg0,'Position');
set(hlg0,'Position',[lgPos(1)+0.05 lgPos(2)-0.225 lgPos(3) lgPos(4)]);
%hyl=ylabel('Depth (mm)')

%hxl=xlabel('Radius');
% xlPos=get(hxl,'Position');
% set(hxl,'Position',[0.9 9.5 xlPos(3)]);

% hyl=ylabel('Depth');
% ylPos=get(hyl,'Position');
% set(hyl,'Position',[0 1.5 ylPos(3)]);
%set(hyl,'Position',[hxl(1) hxl(2) hxl(3) hxl(4)]);
%htext=text(0.025,1,'Depth','FontSize',14);
%hline=line([0 0],[0 9.45]);
%set(hline,'Color','k');
%hline=textarrow('arrow',[0 1],[9.45 9.45]);
%set(hline,'Color','k');
if 1
    
    axes(ha(2));
    set(ha(2),'XAxisLocation','top');
    load('../data/precomputed/RZAfull_092418.mat','R','Z','A');
    ra=R(:); za=Z(:);
    [RR,ZZ] = meshgrid(ra, za);
    hold on
    hcnt1=contourf(RR,ZZ,log10(A),log10(contourLevels),'ShowText','off');
    xl=xlim;
%     SCALP=9.45;
%     SKULL=8.8;
%     CSF=8.2;
%     BRAIN=8;
    plot([xl(1) xl(end)],10*(SCALP-SKULL)*[1 1],'--','Color',get(har4,'FaceColor'))
    plot([xl(1) xl(end)],10*(SCALP-CSF)*[1 1],'--','Color',get(har3,'FaceColor'))
    plot([xl(1) xl(end)],10*(SCALP-BRAIN)*[1 1],'--','Color',get(har2,'FaceColor'))
    %view([90 -90])
    xlabel('Radius (mm)','FontSize',12);
    ylabel('Depth (mm)','FontSize',12);
    caxis(log10([contourLevels(1) contourLevels(length(contourLevels))]));
    colormap(cm)
    %xlim([0 ra(end)]);
    xlim([0 40]);
    ylim([0 40]);
    set(ha(2),'XTick',0:20:40);
    set(ha(2),'Ytick',0:20:40);
    set(ha(2),'XTickLabel',{'0','20','40'});
    set(ha(2),'YtickLabel',{'0','20','40'});
    axis ij;
    axis equal
    x = [0 eps eps 0];
    y = [0 0 eps eps];
    cmInds=[1 33 66 100];
    hpatch1=patch(x,y,cm(cmInds(4),:));
    hpatch2=patch(x,y,cm(cmInds(3),:));
    hpatch3=patch(x,y,cm(cmInds(2),:));
    hpatch4=patch(x,y,cm(cmInds(1),:));
    hlg=legend([hpatch1 hpatch2 hpatch3 hpatch4],{'>10','>1','>0.1','>0.01 J/cm^3'});
    set(hlg,'location','southeast');
    set(hlg,'box','off');
    lgPos=get(hlg,'Position');
    % set(hlg,'FontSize',8);
    set(hlg,'Position',[lgPos(1)+0.02 lgPos(2)+0.1 lgPos(3) lgPos(4)]);
    %set(hlg,'Orientation','horizontal');
    axPos=get(ha(2),'Position');
    set(ha(2),'Position',[axPos(1) axPos(2)+0.1 axPos(3) axPos(4)/1.5]);
    freezeColors
    
    
    
    axes(ha(3))
    %subplot(223);
    hold on
    xmri=-(mrdims(1)-1)/2:(mrdims(1)-1)/2;
    ymri=-(mrdims(2)-1)/2:(mrdims(2)-1)/2;
    [XMRI,YMRI]=ndgrid(1:mrdims(1),1:mrdims(2));
    hpc=pcolor(XMRI,YMRI,anat(:,:,slice2show));
    set(hpc,'EdgeColor','none');
    set(hpc,'FaceAlpha',0.5);
    axis ij;
    axis equal;
    %axis off
    colormap bone;
    freezeColors;
    %logAbsorption(logAbsorption(:,:,slice2show)<0.001)=NaN;
    hpc=pcolor(XMRI,YMRI,logAbsorption(:,:,slice2show));
    set(hpc,'EdgeColor','none');
    set(hpc,'FaceAlpha',0.8);
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
    ylim([10 180]);
    hxl=title(ha(3),'z=88');set(hxl,'Visible','on','FontWeight','normal')
    axis off
    axPos=get(ha(3),'Position');
    set(ha(3),'Position',[axPos(1) axPos(2)+0.05 axPos(3) axPos(4)]);
    
    %%
    % compute total energy absorbed in GM/WM/CSF
%     voxelVolume=0.1^3; % cm^3
%     csfTotal=sum(absorption(labels==1))*voxelVolume;
%     gmTotal=sum(absorption(labels==2))*voxelVolume;
%     wmTotal=sum(absorption(labels==3))*voxelVolume;
    
    axis(ha(4));
    hold(ha(4),'on')
    %hbar2=bar(ha(4),[1 2 3],[csfTotal gmTotal wmTotal],'EdgeColor','none');
    hbar2=bar(ha(4),[1 2 3],[mean(csfTotals) mean(gmTotals) mean(wmTotals)],'EdgeColor','none');    
    hold all
    %xlim([0.75 3.25])
%     plot(ha(4),[1 1],[muCsf-semCsf muCsf+semCsf],'k');
%     plot(ha(4),[2 2],[muGm-semGm muGm+semGm],'k');
%     plot(ha(4),[3 3],[muWm-semWm muWm+semWm],'k');
    %herr=errorbar(ha(4),[muCsf muGm muWm],[semCsf semGm semWm],'LineStyle','none','Color',[0 0 0]);
    herr=errorbar(ha(4),[muCsf muGm muWm],[semCsf semGm semWm]*sqrt(20),'LineStyle','none','Color',[0 0 0]); % std
    set(herr,'CapSize',6)
    set(herr,'LineWidth',2);
    ylabel(ha(4),'Absorbed Energy (J)');
    set(ha(4),'Xtick',1:3);
    set(ha(4),'XTickLabel',{'CSF','GM','WM'});
    set(ha(4),'Ytick',0:0.5:1.5);
    set(ha(4),'YtickLabel',{'0','0.5','1','1.5'});
    axPos=get(ha(4),'Position');
    set(ha(4),'Position',[axPos(1)+0.05 axPos(2)+0.05 axPos(3) axPos(4)]);
    %set(ha(4),'Xlim',[0.75 3.25]);
   
    set(ha(4),'Box','off');
    

    


end

%%
shr=[0.75 0.75 0.75 0.75];
for i=1:4
    pos=get(ha(i),'Position');
    set(ha(i),'Position',[pos(1) pos(2) shr(i)*pos(3) shr(i)*pos(4)]);
end

%%
pos=get(ha(1),'Position');
shr=1.1;
set(ha(1),'Position',[0.025 pos(2) shr*pos(3) shr*pos(4)]);

pos=get(ha(2),'Position');
set(ha(2),'Position',[0.35 pos(2)+0.015 pos(3) pos(4)]);

pos=get(ha(3),'Position');
shr=0.8;
set(ha(3),'Position',[0.58 0.6 shr*pos(3) shr*pos(4)]);

pos=get(ha(4),'Position');
shr=0.6;
set(ha(4),'Position',[0.85 0.65 shr*pos(3) shr*pos(4)]);

lgPos=get(hlg0,'Position');
lgShr=0.8;
set(hlg0,'Position',[0.025 0.59 lgShr*lgPos(3) lgShr*lgPos(4)]);

lgPos=get(hlg,'Position');
lgShr=0.5;
set(hlg,'Position',[0.51 0.6 lgShr*lgPos(3) lgShr*lgPos(4)]);

cbPos=get(hcb,'Position');
shr=0.7;
set(hcb,'Position',[0.725 0.67 shr*cbPos(3) shr*cbPos(4)]);

set(ha(4),'YTick',[0 1]);
lPos=get(get(ha(4),'YLabel'),'Position');
set(get(ha(4),'Ylabel'),'Position',[ -0.95  1.0000   -1.0000]);
%%
sublabel(ha,-15,0,'FontWeight','Bold','FontSize',16);
print -dpng -r600 ../figures/absorptionFigure
crop('../figures/absorptionFigure.png');