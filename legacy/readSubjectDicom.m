clearvars; close all; clc

addpath(genpath('..'));
addpath(genpath('/Users/jacek/Documents/MATLAB/NIfTI_20140122'));
src='/Users/jacek/Documents/MATLAB/LLLT/data/S02/DICOM';
niiFolder='/Users/jacek/Documents/MATLAB/LLLT/data/S02/NII';
fmt=0;
varargout = dicm2nii(src, niiFolder, fmt);

% %%
% 
% nii = load_untouch_nii(fullfile(niiFolder,'BOLD_AP.nii'));
% 
% %%
% data=double(nii.img);
% muData=mean(data,4);
% 
% %%
% data2D=reshape(data,[size(data,1)*size(data,2)*size(data,3) size(data,4)]);
% 
% %%
% nKeep=10000;
% mu_data2D=mean(abs(data2D),2);
% [vals,sortind]=sort(mu_data2D,'descend');
% inds2look=sortind(1:nKeep);
% indsNoise=sortind(end:-1:end-nKeep+1);
% 
% peakVoxels=data2D(inds2look,:);
% lowVoxels=data2D(indsNoise,:);
% 
% 
% %%
% figure;
% subplot(211)
% plot(mean(peakVoxels,1))
% hold on
% %plot(mean(lowVoxels,1),'k')
% subplot(212);
% plot(mean(lowVoxels,1),'k')
% %plot(smoothts(mean(peakVoxels,1))); ylim([700 800])
% 
% %%
% figure;






