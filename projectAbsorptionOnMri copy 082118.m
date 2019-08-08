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
load '../data/mcml/forJacek_trimmed';
axialCutOff=axialCutOff*10; % mm
radialCutOff=radialCutOff*10; % mm
xAxis = linspace(0,axialCutOff, size(thinnedTruncatedAverageData,1));
yAxis = linspace(0,radialCutOff, size(thinnedTruncatedAverageData,1));
A=thinnedTruncatedAverageData;
[X,Y] = meshgrid(xAxis, yAxis);
XX=cat(2,-X(:,end:-1:1),X);
YY=[Y,Y];
AA=cat(2,A(:,end:-1:1),A);


%%

% resample on a 1 mm grid instead
% dsr=10; % from 0.1 mm to 1 mm
% [Xq,Yq]=meshgrid(resample(xAxis,1,dsr),resample(xAxis,1,dsr));
% Xq=cat(2,-Xq(:,end:-1:1),Xq);
% Yq=[Yq,Yq];
% Aq = interp2(XX,YY,AA,Xq,Yq);

% % shift to laser center
% Xq=Xq;
% Yq=Yq-85;

theta=atan(9/27); % measured from afni
XXn=XX*cos(theta)+YY*sin(-theta);
YYn=XX*sin(theta)+YY*cos(theta);

XXn=XXn+laserOrigin(1);
YYn=YYn+laserOrigin(2);

XX=XXn;
YY=YYn;

%% draw
figure;
subplot(221); hold on
hcnt1=contourf(XX,YY,log10(AA),log10(contourLevels),'ShowText','of');
%hcnt2=contourf(-XX,Y,log10(A),log10(contourLevels),'ShowText','off');
%xlim([-30 30]);
%ylim([0 30]);
xlabel('Radial distance (cm)');
ylabel('Axial distance (cm)');
caxis(log10([contourLevels(1) contourLevels(length(contourLevels))]));
axis ij; axis equal
colormap(cm)

hs(2)=subplot(222); 
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
%%
figure;
hold on
[C1,hcnt1]=contourf(XX,YY,log10(AA),log10(contourLevels),'ShowText','off');
caxis(log10([contourLevels(1) contourLevels(length(contourLevels))]));
colormap(cm);
axis ij; 
axis equal; 
xmri=-(mrdims(1)-1)/2:(mrdims(1)-1)/2;
ymri=-(mrdims(2)-1)/2:(mrdims(2)-1)/2;
[XMRI,YMRI]=ndgrid(1:mrdims(1),1:mrdims(2));
freezeColors;
hpc=pcolor(XMRI,YMRI,anat(:,:,slice2show));
set(hpc,'EdgeColor','none');
colormap bone;






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

