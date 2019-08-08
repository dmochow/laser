clear variables; close all; clc
% run afni preprocessing from matlab

% specify the subject here
subjStr='S10';
pathToData=['../data/' subjStr '/NII/'];

!export PATH=$PATH:/Users/jacekdmochowski/abin

cd(pathToData);


str1='!@auto_tlrc -base TT_N27+tlrc -input anat.nii';
eval(str1);

%
%str2='!@auto_tlrc -apar anat_at.nii -input roiMask.nii';
str2='!@auto_tlrc -apar anat_at.nii -input roi+orig.nii';
eval(str2);

% 
str3='!align_epi_anat.py -anat anat.nii -epi bold_e2.nii -epi_base 5 -child_epi bold_e1.nii bold_e3.nii -epi2anat -suffix _al_at -tlrc_apar anat_at.nii';
eval(str3);
% % now align the epis to the anatomy

% resample the mask on the bold grid
str4='!3dresample -master bold_e2_tlrc_al_at+tlrc -prefix roiMask_at_resampled -input roiMask_at.nii';
eval(str4);

%resample the atlased anatomy on the bold grid? (try to get the brain mask with
% this)
str5='!3dresample -master bold_e2_tlrc_al_at+tlrc -prefix anat_at_resampled -input anat_at.nii';
eval(str5);

%
% % resample mask on the aligned and atlased bold
% !3dresample -master bold_e2_tlrc_al_at+tlrc -prefix roiMask_at_resampled -input roiMask_at.nii
