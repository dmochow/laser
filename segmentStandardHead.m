clear all; close all; clc;
!export PATH=$PATH:/Users/jacekdmochowski/abin
PATH = getenv('PATH');
setenv('PATH', [PATH ':/Users/jacekdmochowski/abin']);


basePath=['../data/SAVG/NII/'];

% segment skull stripped anatomical
str=['!3dSeg -anat ' basePath 'TT_N27+tlrc -mask AUTO -classes ''CSF ; GM ; WM'' -bias_classes ''GM ; WM'' -bias_fwhm 25 -mixfrac UNI -main_N 5 -blur_meth BFT -prefix ' basePath];
eval(str);