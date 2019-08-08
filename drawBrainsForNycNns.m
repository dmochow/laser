clear all; close all; clc
% draw some images for nyc nans 2018

pathToData='../data/SAVG/NII';
anatFilename=fullfile(pathToData,'TT_N27+tlrc');
brainMaskFilename=fullfile(pathToData,'brain_mask+tlrc.BRIK');
roiMaskFilename=fullfile(pathToData,'mu_roi_r19_z26+tlrc.BRIK');
sliceIndex=92;
cropMri=1;

[~, anat, ~, ~] = BrikLoad (anatFilename);
anat=anat(end:-1:1,:,:);

[~, brainMask, ~, ~] = BrikLoad (brainMaskFilename);
brainMask=brainMask(end:-1:1,:,:);

[~, roiMask, ~, ~] = BrikLoad (roiMaskFilename);
roiMask=roiMask(end:-1:1,:,:);
roiMask=logical(roiMask*0);

% create dummy variables
isSig1=zeros(size(anat))==1;
isSig2=zeros(size(anat))==1;

h = colorBrain(anat(:,:,sliceIndex),roiMask(:,:,sliceIndex),isSig1(:,:,sliceIndex)==1,brainMask(:,:,sliceIndex),cropMri);
print -dpng ../figures/anatBrain