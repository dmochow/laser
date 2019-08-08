clear all; close all; clc

%% coordinates of average laser origin
% extracted by visual inspection of muroi
% in Afni:
% sagittal: 104
% coronal: 5 (225-5)
% axial: 87
% angle atan(9/27);

contourLevels=[0.01 0.1 1 10];

%% colormap
cmap = cbrewer('div','PuOr',100);
cm=flipud(cmap);

%% anatomy
pathToData='../data/SAVG/NII';
anatFilename=fullfile(pathToData,'TT_N27+tlrc');
brainMaskFilename=fullfile(pathToData,'brain_mask+tlrc.BRIK');
[~, anat, info, ~] = BrikLoad (anatFilename);
anat=anat(end:-1:1,:,:);
[~, brainMask, ~, ~] = BrikLoad (brainMaskFilename);
brainMask=brainMask(end:-1:1,:,:);
slice2show=88;
laserOrigin=[104 5];
mrdims=size(anat);
anat(~brainMask)=NaN;

%% absorption
vars=load('../data/mcml/forJacek_trimmed');
axialCutOff=axialCutOff*10; % mm
radialCutOff=radialCutOff*10; % mm
xAxis = linspace(0,axialCutOff, size(thinnedTruncatedAverageData,1));
yAxis = linspace(0,radialCutOff, size(thinnedTruncatedAverageData,1));
A=thinnedTruncatedAverageData;
[X,Y] = meshgrid(xAxis, yAxis);
XX=cat(2,-X(:,end:-1:1),X);
YY=[Y,Y];
AA=cat(2,A(:,end:-1:1),A);
XXo=XX; YYo=YY; % original values (unrotated, unshifted)

%% write raw absorption to disk
R=X; Z=Y; 

% extend to 3D
nPages=10;
R3=repmat(R,[1 1 nPages]);
Z3=repmat(Z,[1 1 nPages]);
A3=repmat(A,[1 1 nPages]);
TH3=[];
TH3(1,1,:)=linspace(0,2*pi,nPages);
TH3=repmat(TH3,[300 300 1]);
[X3,Y3,Z3] = pol2cart(TH3,R3,Z3); 


% now interpolate
F = scatteredInterpolant(X3(:),Y3(:),Z3(:),A3(:));
[xq,yq,zq]=ndgrid(linspace(-max(X3(:)),max(X3(:)),21), linspace(-max(Y3(:)),max(Y3(:)),21) , linspace(0,max(Z3(:)),10) );
Aq = F(xq,yq,zq);
%%
figure
for sl=1:10
    subplot(2,5,sl);
    imagesc(Aq(:,:,sl)); axis square
end

%%

figure
for sl=1:10
    subplot(2,5,sl);
    imagesc(squeeze(Aq(:,sl,:))); axis square
end

figure
for sl=1:10
    subplot(2,5,sl);
    imagesc(squeeze(Aq(sl,:,:))); axis square
end


%%




%%
% nSlicesToAdd=300;
% A3=zeros(size(AA,2),1,size(AA,1)); % r, theta, z
% A3(:,1,:)=permute(AA,[2 1]);
% A3=repmat(A3,[1 nSlicesToAdd 1]);
% 
% figure; 
% for s=1:5
%     subplot(3,5,s); imagesc(squeeze(A3(:,10*s,:))); axis square
%     subplot(3,5,s+5); imagesc(squeeze(A3(:,:,10*s))); axis square
%     subplot(3,5,s+10); imagesc(squeeze(A3(10*s,:,:))); axis square
% end
% 

%save('../data/precomputed/RZA.mat','R','Z','A');
%%


theta=atan(9/27); % measured from afni
XXn=XX*cos(theta)+YY*sin(-theta);
YYn=XX*sin(theta)+YY*cos(theta);

XXn=XXn+laserOrigin(1);
YYn=YYn+laserOrigin(2);

XX=XXn;
YY=YYn;

%% draw
figure;
%subplot(221); 
hold on
hcnt1=contourf(XXo,YYo,log10(AA),log10(contourLevels),'ShowText','of');
xlabel('Lateral extent (mm)','FontSize',20);
ylabel('Depth (mm)','FontSize',20);
caxis(log10([contourLevels(1) contourLevels(length(contourLevels))]));
%axis ij; 
axis equal
colormap(cm)
x = [0 eps eps 0];
y = [0 0 eps eps];
cmInds=[1 33 66 100];
hpatch1=patch(x,y,cm(cmInds(4),:));
hpatch2=patch(x,y,cm(cmInds(3),:));
hpatch3=patch(x,y,cm(cmInds(2),:));
hpatch4=patch(x,y,cm(cmInds(1),:));

% add sidebar showing layers
SCALP=6;
SKULL=11;
CSF=11.1;
BRAIN=30;

%area([30 40],[SKULL SKULL],'BaseValue',SCALP);

area([30 40],[BRAIN BRAIN],'BaseValue',CSF,'FaceAlpha',0.5,'LineStyle','none');
area([30 40],[CSF CSF],'BaseValue',SKULL,'FaceAlpha',0.5,'LineStyle','none');
area([30 40],[SKULL SKULL],'BaseValue',SCALP,'FaceAlpha',0.5,'LineStyle','none');
area([30 40],[SCALP SCALP],'BaseValue',0,'FaceAlpha',0.5,'LineStyle','none'); 
text(32,SCALP/2,'SCALP');
text(32,SCALP+(SKULL-SCALP)/2,'SKULL');
text(32,SKULL+(BRAIN-CSF)/2,'BRAIN');
text(32,SKULL+(BRAIN-CSF)/2,'BRAIN');
annotation('textarrow',[0.92 0.9],[0.46 0.46],'String','CSF')

hlg=legend([hpatch1 hpatch2 hpatch3 hpatch4],{'10 J','1 J','0.1 J','0.01 J'},'FontSize',16);
lgPos=get(hlg,'Position');
set(hlg,'Position',[lgPos(1)-0.3 lgPos(2)+0.175 lgPos(3) lgPos(4)]);
set(hlg,'Orientation','horizontal');

set(gca,'Xtick',[-30:5:30]);

print -dpng -r600 ../figures/figure0
crop('../figures/figure0.png');



%%
figure
hold on
xmri=-(mrdims(1)-1)/2:(mrdims(1)-1)/2;
ymri=-(mrdims(2)-1)/2:(mrdims(2)-1)/2;
[XMRI,YMRI]=ndgrid(1:mrdims(1),1:mrdims(2));
hpc=pcolor(XMRI,YMRI,anat(:,:,slice2show));
set(hpc,'EdgeColor','none');
axis ij; 
axis equal; 
colormap bone;
freezeColors;
[C1,hcnt1]=contourf(XX,YY,log10(AA),log10(contourLevels),'ShowText','off');
caxis(log10([contourLevels(1) contourLevels(length(contourLevels))]));
axis equal
axis off
colormap(cm);
freezeColors

print -dpng ../figures/figure1
crop('../figures/figure1.png');
% %%
% figure;
% hold on
% [C1,hcnt1]=contourf(XX,YY,log10(AA),log10(contourLevels),'ShowText','off');
% caxis(log10([contourLevels(1) contourLevels(length(contourLevels))]));
% colormap(cm);
% axis ij; 
% axis equal; 
% xmri=-(mrdims(1)-1)/2:(mrdims(1)-1)/2;
% ymri=-(mrdims(2)-1)/2:(mrdims(2)-1)/2;
% [XMRI,YMRI]=ndgrid(1:mrdims(1),1:mrdims(2));
% freezeColors;
% hpc=pcolor(XMRI,YMRI,anat(:,:,slice2show));
% set(hpc,'EdgeColor','none');
% colormap bone;
% 
% %%
% figure





%% JUNKHEAP
% subplot(222); hold on
% hcnt1=contourf(Xq,Yq,log10(Aq),log10(contourLevels),'ShowText','off');
% hcnt2=contourf(-Xq,Yq,log10(Aq),log10(contourLevels),'ShowText','off');
% xlabel('Radial distance (cm)');
% ylabel('Axial distance (cm)');
% caxis(log10([contourLevels(1) contourLevels(length(contourLevels))]));
% axis equal


% fillColors=[1 0 0; 0 1 0; 0 0 1; 1 1 0];
% s1=contourdata(C1);
% s2=contourdata(C2);
% hs(4)=subplot(224); hold on
% for i=1:numel(s1)
%     
%     [~,ind1]=min(abs(s1(i).level-log10(contourLevels)));
%     fill(s1(i).xdata,s1(i).ydata,fillColors(ind1,:),'FaceAlpha',0.5);
%     
% %     [~,ind2]=min(abs(s2(i).level-log10(contourLevels)));
% %     fill(s2(i).xdata,s2(i).ydata,fillColors(ind2,:),'FaceAlpha',0.5);
%    
% end
% caxis(log10([contourLevels(1) contourLevels(length(contourLevels))]));
% axis equal

% %[C,h] = contourf(X,Y,log10(A),[log10(contourLevels(4)) log10(contourLevels(4)) ]); 
% [C,h] = contourf(log10(A),'LevelList',log10(contourLevels(1))); 
% % Close the contour by adding two extra points
% % x = [C(1,2:end), X(end), X(end)];
% % y = [C(2,2:end), Y(1), Y(end)];
% x = [C(1,:)];
% y = [C(2,:)];
% fill(x,y,'m','FaceAlpha',0.5);

