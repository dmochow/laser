clear all; close all; clc

pathToData='../data/thermo/scanner_data/tS01a';
niiFilename='anat.nii';
nii=load_untouch_nii(fullfile(pathToData,niiFilename));

img=nii.img;
img2=double(nii.img);

%%
save(fullfile(pathToData,'anat_am.mat'),'img','img2');

%%
load(fullfile(pathToData,'anat_am.mat'),'img','img2');

%%
size(img2)
figure;
subplot(2,2,1);
imagesc(squeeze(img2(:,:,180))); colormap bone
subplot(2,2,2);
imagesc(squeeze(img2(:,120,:))); colormap bone
subplot(2,2,3);
imagesc(squeeze(img2(104,:,:)).'); axis xy; colormap bone