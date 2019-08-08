clear all; close all; clc

subjNum=7;
afniCompIndx=1;


betasFilename=['/Users/jacek/Documents/MATLAB/LLLT/data/S0' num2str(subjNum) '/NII/mebold2go/meica.bold_e123/TED/betas_OC.nii'];
mixFilename=['/Users/jacek/Documents/MATLAB/LLLT/data/S0' num2str(subjNum) '/NII/mebold2go/meica.bold_e123/TED/meica_mix.1D'];

nii=load_untouch_nii(betasFilename);
betas=nii.img;
sliceIndx=[41,7,24]+[1 1 1];

A=importdata(mixFilename,' ');


figure;
subplot(221)
imagesc(squeeze(betas(sliceIndx(1),:,:,afniCompIndx+1))); colormap bone
%subplot(222)
%imagesc(squeeze(betas(:,sliceIndx(2),:,afniCompIndx+1))); colormap bone
%subplot(223)
%imagesc(squeeze(betas(:,:,sliceIndx(3),afniCompIndx+1))); colormap bone


subplot(224);
plot(A(:,afniCompIndx+1),'k'); % afni starts at 0, matlab at 1

%figFilename=['../figures/s05_mixica_' num2str(afniCompIndx)];
%print('-depsc',figFilename);