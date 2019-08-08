clear all; close all; clc

%% CONSTANTS
a=-0.01/1e6; % "a" Yuan et al. 2012
gamma=42.58e6; % gyromagnetic ratio of hydrogen in water
Bo=3; % 3 Tesla
TE=17e-3; % 17 ms?
delrad=8*pi/(2^16-1); % 
gmIndex=2; 

% 
subjStr='tS01';
pathToDicom=['../data/thermo/' subjStr '/'];
pathToData=['../data/thermo/output/' subjStr '/'];
pathToRefData=['../data/thermo/output/' subjStr '/'];
pathToMagData=['../data/thermo/output/' subjStr '/'];

% get dicom info
files=dir(fullfile(pathToDicom,'*.ima'));
info = dicominfo(fullfile(pathToDicom,files(1).name));
slope=info.RescaleSlope;
intercept=info.RescaleIntercept;

inputPhFilename=fullfile(pathToData,['dsph_tlrc_al+tlrc.BRIK']);
[~, inputPh, iInfo, ~] = BrikLoad (inputPhFilename);
temp=inputPh/(a*gamma*Bo*TE);

roiMaskFilename=fullfile(pathToData,['resampled_roi_r21_z39+tlrc.BRIK']);
[~, roiMask, ~ , ~] = BrikLoad (roiMaskFilename);

% grey matter mask
resampledSegFilename=fullfile(pathToData,'resampled_Classes+tlrc.BRIK');
[~, seg, ~, ~] = BrikLoad (resampledSegFilename);
greyMask=logical(seg==gmIndex);

finalMask=roiMask&greyMask;
roiTemp=vol2ts(temp,finalMask);

%%
regMask=~roiMask&greyMask;
regTemp=vol2ts(temp,regMask);
[U,S,V]=svd(regTemp,0);
tempClean=regressOut(roiTemp.',U(:,1:10).',1).';

%x=double(intercept+double(x)*slope);



%%
% convert to temperature change


