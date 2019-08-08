clear variables; close all; clc
% run afni preprocessing from matlab

% specify the subject here
subjStr='S04';
pathToData=['../data/' subjStr '/NII/'];

%% FINAL DATASETS
% resampled_anat+tlrc
% resampled_roi+tlrc
% bold_ex_tlrc_al+tlrc


% do the following in MATLAB on a new machine
PATH = getenv('PATH');
setenv('PATH', [PATH ':/Users/jacekdmochowski/abin']);

!export PATH=$PATH:/Users/jacekdmochowski/abin

cd(pathToData);

%% PIPELINE 1: with atlasing to talairach
str0='!3dcopy anat.nii anat';
eval(str0);

str0='!3dcopy bold_e1.nii bold_e1';
eval(str0);

str0='!3dcopy bold_e2.nii bold_e2';
eval(str0);

str0='!3dcopy bold_e3.nii bold_e3';
eval(str0);

str1='!@auto_tlrc -base TT_N27+tlrc -input anat+orig';
eval(str1);
% align the anatomical to the talairach atlas

%
%str2='!@auto_tlrc -apar anat_at.nii -input roiMask.nii';
str2='!@auto_tlrc -apar anat+tlrc -input roi+orig';
eval(str2);
% align the roi mask to the talairach atlas

% 
str3='!align_epi_anat.py -anat anat+orig -epi bold_e2+orig -epi_base 5 -child_epi bold_e1+orig bold_e3+orig -epi2anat -suffix _al -tlrc_apar anat+tlrc';
eval(str3);
% % now align the epis to the anatomy, and also telling it that we have an
% atlased anatomy

% resample the mask on the bold grid
str4='!3dresample -master bold_e2_tlrc_al+tlrc -prefix resampled_roi -input roi+tlrc';
eval(str4);
% 
% %resample the atlased anatomy on the bold grid? (try to get the brain mask with
% % this)
str5='!3dresample -master bold_e2_tlrc_al+tlrc -prefix resampled_anat -input anat+tlrc';
eval(str5);
% % 

