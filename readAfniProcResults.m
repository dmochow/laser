clear all; close all; clc
addpath(genpath('/Users/jacekdmochowski/PROJECTS/COMMON'));
!export PATH=$PATH:/Users/jacekdmochowski/abin
PATH = getenv('PATH');
setenv('PATH', [PATH ':/Users/jacekdmochowski/abin']);
subjNum=23;
subjStr=['S' num2str(subjNum)];
pathToData=['../data/' subjStr '/' subjStr '.results'];

epiFilename=['all_runs.' subjStr '+tlrc.BRIK'];
[~, bold, bInfo, ~] = BrikLoad (fullfile(pathToData,epiFilename));

maskFilename=['mask_epi_anat.' subjStr '+tlrc.BRIK'];
[~, mask, mInfo, ~] = BrikLoad (fullfile(pathToData,maskFilename));

laserFilename=['follow_ROI_LASER+tlrc.BRIK'];
[~, laser, lInfo, ~] = BrikLoad (fullfile(pathToData,laserFilename));

% %%
% origpath=pwd;
% cd(pathToData);
% str=['!3dresample -master all_runs.' subjStr '+tlrc -prefix resampled_Classes -input Classes+tlrc']; eval(str);
% cd(origpath);
% %%
% 
% segFilename=['resampled_Classes+tlrc.BRIK'];
% [~, seg, sInfo, ~] = BrikLoad (fullfile(pathToData,segFilename));

finalMask=laser==1 & mask>0;
%finalMask=laser==1 & seg==2;
roits=vol2ts(bold,finalMask);

figure;
plot(mean(roits,2));
