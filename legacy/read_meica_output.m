clear all; close all; clc

% read in the output of MEICA

% point the following line of code to where you unzipped NIfTI matlab
% package
addpath(genpath('/Users/jacek/Documents/MATLAB/NIfTI_20140122'));

% point the program to where the denoise data lies
dataFilename=['../data/S05/NII/mebold2go/bold_e123_medn.nii'];

% point the program to where the brain mask lies (so that we only look at
% brain voxels)
maskFilename=['../data/S05/NII/mebold2go/meica.bold_e123/eBvrmask.nii'];

% read in mask
maskNii=load_untouch_nii(maskFilename);
mask=logical(maskNii.img);
figure;imagesc(mask(:,:,30)); colormap bone

% read in data and only look at brain
dataNii=load_untouch_nii(dataFilename);
data=double(dataNii.img);
data_4D=data;
data=reshape(data,[size(data,1)*size(data,2)*size(data,3) size(data,4)]);
brainData=data(mask(:),:);

% show the time course of a single voxel
voxelIndex=50000;
figure
plot(brainData(voxelIndex,:));











