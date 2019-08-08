clear all; close all; clc

% anatomy
pathToData='../data/SAVG/NII';
anatFilename=fullfile(pathToData,'TT_N27+tlrc');
[~, anat, info, ~] = BrikLoad (anatFilename);
anat=anat(end:-1:1,:,:);
slice2show=100;
laserOrigin=[20 40];

% absorption
%load('../data/mcml/Arz','Abs','r','z');
%A=Abs.'; clear Abs;
load '../data/mcml/forJacek_trimmed';
axialCutOff=axialCutOff*10; % mm
radialCutOff=radialCutOff*10; % mm
xAxis = linspace(0,axialCutOff, size(thinnedTruncatedAverageData,1));
yAxis = linspace(0,radialCutOff, size(thinnedTruncatedAverageData,1));
hLines = [0.6, 1.1, 1.2];
[X,Y] = meshgrid(xAxis, yAxis);
A=thinnedTruncatedAverageData;

% resample on a 1 mm grid instead
dsr=10; % from 0.1 mm to 1 mm
%[Z,R]=meshgrid(r,z);
mrdims=[161,191,151];
[Xq,Yq]=meshgrid(0:1:80,0:1:190);
Aq = interp2(X,Y,A,Xq,Yq);
%Aqs=circshift(Aq,[laserOrigin(1) laserOrigin(2)]);
%%
figure;
subplot(221); hold on
hcnt1=contour(X,Y,A,[0.1 1 10 100]);
hcnt2=contour(-X,Y,A,[0.1 1 10 100]);
xlim([-30 30]);
ylim([0 30]);
xlabel('Radial distance (cm)');
ylabel('Axial distance (cm)');

subplot(222);
imagesc(anat(:,:,slice2show));

% subplot(223);
% hpc2=pcolor(Rq,Zq,Aq);
% shading interp
% 
% subplot(224);
% hpc2=pcolor(Rq,Zq,Aqs);
% shading interp



% figure;
% subplot(221);
% imagesc(Abs);
% subplot(222);
% imagesc(y,x,Abs);xlabel('depth'); ylabel('radius (lateral)');
% subplot(223);
% contour(z,r,Abs,[0.001 0.001]);

