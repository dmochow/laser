clearvars; close all; clc

addpath(genpath('..'));

niiFilename='/Users/jacek/Documents/MATLAB/LLLT/data/S01/NII/suBOLD_AP.nii';
nii = load_untouch_nii(niiFilename);
bold=double(nii.img);
bold_2D=reshape(bold,[size(bold,1)*size(bold,2)*size(bold,3) size(bold,4)]);

%% simple ratio of powers before vs during laser
pwr_pre=mean(bold_2D(:,1:400),2);
pwr_laser=mean(bold_2D(:,401:800),2);
pwr_ratio=(pwr_laser-pwr_pre)./pwr_pre;

%%
pwr_ratio_3D=reshape(pwr_ratio,[size(bold,1) size(bold,2) size(bold,3)]);

%%
figure
subplot(211);
vol3d('cdata',bold(:,:,:,1));
subplot(212);
h = vol3d('cdata',pwr_ratio_3D);

%%
