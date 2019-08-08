clear all; close all; clc
addpath(genpath('/Users/jacekdmochowski/PROJECTS/COMMON'));

%% PARAMETERS
subjStr='tS01c';
nTRs=150;
TR=20.5*60/nTRs; % HARD-CODED FOR NOW
onsetTR=21; % 
offsetTR=21+round(10*60/TR);
P=4; Q=4; zero=1; override=1; % params for bad voxels
nPCs=3;
% 
% nTRs=50;
% TR=20.5*60/nTRs; % HARD-CODED FOR NOW
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


% 
% regress out motion
% create the motion regressors
motionFilename=fullfile(pathToData,['mag_vr_motion.1D']);
fid=fopen(motionFilename);
mals=textscan(fid,'%f %f %f %f %f %f');
mals=cell2mat(mals);
dmals=cat(1,zeros(1,6),diff(mals));
mals=cat(2,mals,dmals);
mals=mals-repmat(mean(mals),size(mals,1),1);
iwBoldTs=vol2ts(tempChange);
owBoldTs = regressOut(iwBoldTs.',mals.',1).';
tempChange=ts2vol(owBoldTs,dummyMask);

roiMaskFilename=fullfile(pathToData,['resampled_roi_r21_z39+orig.BRIK']);
[~, roiMask, ~ , ~] = BrikLoad (roiMaskFilename);

% grey matter mask
resampledSegFilename=fullfile(pathToData,'resampled_Classes+orig.BRIK');
[~, seg, ~, ~] = BrikLoad (resampledSegFilename);
greyMask=logical(seg==gmIndex);

finalMask=roiMask&greyMask;
roiTempChange=vol2ts(tempChange,finalMask);

%bring in drift file
fid=fopen(fullfile(pathToScannerData,'drift.txt'));
formatSpec=['rep %u phase = %f \n'];
fgets(fid); fgets(fid);
D = fscanf(fid,formatSpec);
fclose(fid);
drift=D(2:2:end);

% regress out drift
roiTempChange=regressOut(roiTempChange.',drift',1).';

%%
regMask=~roiMask&greyMask;
regTemp=vol2ts(tempChange,regMask);
[U,S,V]=svd(regTemp,0);
roiTempClean=regressOut(roiTempChange.',U(:,1:nPCs).',1).';

% % remove bad voxels
% roiTempClean=nanBadChannels(roiTempClean.',P,Q,zero,override).';
% 
% % remove bad samples
% roiTempClean = nanBadSamples(roiTempClean.',P,Q,zero,1,1).';


%%
time=(0:nTRs-1)*TR;
figure; hold on
plot(time,mean(roiTempChange,2),'LineWidth',1.5);
plot(time,mean(roiTempClean,2),'LineWidth',1.5);
yl=ylim;
area([onsetTR*TR offsetTR*TR],[yl(2) yl(2)],'BaseValue',yl(1),'FaceColor',[0.6 0.6 0.6],'FaceAlpha',0.15,'LineStyle','none');
text(onsetTR*TR+210,yl(1)+1,'Laser on','FontSize',16);
legend('Raw','After PCA regression');
ylabel('Temp Change (^oC)');
xlabel('Time (s)');
axis tight
print('-dpng',['../figures/thermo_' subjStr]);
crop(['../figures/thermo_' subjStr '.png'],0);
%x=double(intercept+double(x)*slope);



%%
% convert to temperature change


