clear all; close all; clc

% read in the output of MEICA

% point the following line of code to where you unzipped NIfTI matlab
% package
addpath(genpath('/Users/jacek/Documents/MATLAB/NIfTI_20140122'));

% point the program to where the denoise data lies
dataFilename=['../data/S05/NII/mebold2go/anat.nii'];

% % point the program to where the brain mask lies (so that we only look at
% % brain voxels)
% maskFilename=['../data/S05/NII/mebold2go/meica.bold_e123/eBvrmask.nii'];


%%

% read in data from volume
dataNii=load_untouch_nii(dataFilename);
data=dataNii.img;
figure;
subplot(221);
imagesc(squeeze(data(:,:,120))); colormap bone
view([-90 90])
subplot(222);
imagesc(squeeze(data(:,120,:))); colormap bone
view([-90 90])
subplot(223);
imagesc(squeeze(data(120,:,:))); colormap bone
view([-90 90])


%%
figure;
vol3d('Cdata',data);






