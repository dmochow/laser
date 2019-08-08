clear all; close all; clc
% look at the spectrum of bold and try to design the right lowpass filter
addpath(genpath('/Users/jacek/Documents/MATLAB/NIfTI_20140122'));
dataFilename='../data/S01/FSL.orig/BOLD_AP.ica/filtered_func_data.nii';
maskFilename='../data/S01/FSL.orig/BOLD_AP.ica/filtered_func_data.ica/mask.nii';
fs=1/1.5; % 1/TR
maskNii=load_untouch_nii(maskFilename);
mask=logical(maskNii.img);
dataNii=load_untouch_nii(dataFilename);
data=double(dataNii.img);
data=reshape(data,[size(data,1)*size(data,2)*size(data,3) size(data,4)]);
brainData=data(mask(:),:);

%%
nfft=2^nextpow2(size(brainData,2));
brainDataFft=fft(brainData,nfft,2);

%%
freqs=(0:nfft-1)/nfft*fs;
powerSpectrum=brainDataFft.*conj(brainDataFft);
muPowerSpectrum=mean(powerSpectrum,1);
figure;
plot(freqs,muPowerSpectrum);


