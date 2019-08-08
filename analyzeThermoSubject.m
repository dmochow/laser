clear all; close all; clc
addpath(genpath('/Users/jacekdmochowski/PROJECTS/COMMON'));

%% PARAMETERS
subjStr='tS03b';
nPCs=3;

nTRs=150;
TR=20.5*60/nTRs; % HARD-CODED FOR NOW
onsetTR=21; % 
offsetTR=21+round(10*60/TR);

% 
% nTRs=50;
% TR=10.25*60/nTRs; % HARD-CODED FOR NOW
% onsetTR=6; % 
% offsetTR=43;
%%

%% CONSTANTS
a=-0.01e-6; % "a" Yuan et al. 2012
gamma=2*pi*42.58e6; % gyromagnetic ratio of hydrogen in water
Bo=3; % 3 Tesla
TE=17e-3; % 17 ms?
delrad=8*pi/(2^16-1); % 
gmIndex=2; 
nX=64; nY=64; nZ=32; % dimensions of mag/ph
dummyMask=true(nX,nY,nZ); % used for going from 2D to 4D
slope=8*pi/2^16;
intercept=-4*pi;

pathToData=['../data/thermo/output/' subjStr '/'];
pathToScannerData=['../data/thermo/scanner_data/' subjStr '/'];

inputPhFilename=fullfile(pathToData,['smph_al+orig.BRIK']);
[~, inputPh, iInfo, ~] = BrikLoad (inputPhFilename);
inputPh=inputPh/4;  % 08.05.19: compress to (-pi,pi) range from (-4pi,4pi)
tempChange=inputPh/(a*gamma*Bo*TE); % PRF formula

%bring in drift file
fid=fopen(fullfile(pathToScannerData,'drift.txt'));
%fid=fopen(fullfile(pathToScannerData,'drift2.txt'));
formatSpec=['rep %u phase = %f \n'];
fgets(fid); fgets(fid);
D = fscanf(fid,formatSpec);
fclose(fid);
drift=D(2:2:end);

% mask
roiMaskFilename=fullfile(pathToData,['resampled_roi_r21_z39+orig.BRIK']);
[~, roiMask, ~ , ~] = BrikLoad (roiMaskFilename);

resampledSegFilename=fullfile(pathToData,'resampled_Classes+orig.BRIK');
[~, seg, ~, ~] = BrikLoad (resampledSegFilename);
greyMask=logical(seg==gmIndex);

% grab time series of roi voxels
finalMask=roiMask&greyMask;
roiTempChange=vol2ts(tempChange,finalMask);

% regress out drift
[r,p]=corrcoef(mean(roiTempChange,2),drift)
roiTempChange=regressOut(roiTempChange.',drift',1).';


% % regress out motion
%motionFilename=fullfile(pathToData,['mag_vr_motion.1D']);
motionFilename=fullfile(pathToData,['ph_vr_motion.1D']);
fid=fopen(motionFilename);
mals=textscan(fid,'%f %f %f %f %f %f');
mals=cell2mat(mals);
dmals=cat(1,zeros(1,6),diff(mals));
mals=cat(2,mals,dmals);
mals=mals-repmat(mean(mals),size(mals,1),1);
roiTempChange = regressOut(roiTempChange.',mals.',1).';


%%
regMask=~roiMask&greyMask;
regTemp=vol2ts(tempChange,regMask);
[U,S,V]=svd(regTemp,0);
roiTempClean=regressOut(roiTempChange.',U(:,1:nPCs).',1).';

%%
time=(0:nTRs-1)*TR;
figure; hold on
plot(time,mean(roiTempChange,2),'LineWidth',1.5);
plot(time,mean(roiTempClean,2),'LineWidth',1.5);
ylim([-1 1]);
yl=ylim;
area([onsetTR*TR offsetTR*TR],[yl(2) yl(2)],'BaseValue',yl(1),'FaceColor',[0.6 0.6 0.6],'FaceAlpha',0.15,'LineStyle','none');
%text(onsetTR*TR+210,yl(1)+1,'Laser on','FontSize',16);
legend('Raw','After PCA regression');
ylabel('Temp Change (^oC)');
xlabel('Time (s)');
axis tight
print('-dpng',['../figures/thermo_' subjStr]);
crop(['../figures/thermo_' subjStr '.png'],0);


