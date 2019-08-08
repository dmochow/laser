clearvars; close all; clc

addpath(genpath('..'));

niiFilename='/Users/jacek/Documents/MATLAB/LLLT/data/S01/NII/T1w_MPR_BIC_v1.nii';
nii = load_untouch_nii(niiFilename);
data_t1=double(nii.img);

niiFilename='/Users/jacek/Documents/MATLAB/LLLT/data/S01/NII/rBOLD_AP.nii';
nii = load_untouch_nii(niiFilename);
data_bold=double(nii.img);


%%
firstVolume = squeeze(data(:,:,:,1));
h = vol3d('cdata',firstVolume);


guessIndx=[70 35 40];
timeSeries=squeeze(data(guessIndx(1),guessIndx(2),guessIndx(3),:));
figure;
subplot(211);
plot(timeSeries);


guessIndx=[30 45 30];
subplot(212);
timeSeries=squeeze(data(guessIndx(1),guessIndx(2),guessIndx(3),:));
plot(timeSeries);