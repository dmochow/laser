clear all; close all; clc;
!export PATH=$PATH:/Users/jacekdmochowski/abin
PATH = getenv('PATH');
setenv('PATH', [PATH ':/Users/jacekdmochowski/abin']);
subjStrs={'S02','S03','S04','S05','S06','S07','S09','S10','S11','S12','S13','S14',...
    'S15','S16','S17','S18','S19','S20','S21','S22'};
nSubjects=numel(subjStrs);

for s=1:nSubjects
    subjIndx=subjStrs{s};
    basePath=['../data/' subjIndx '/NII/'];
    
    % segment skull stripped anatomical
    str=['!3dSeg -anat ' basePath 'anat_ns+orig -mask AUTO -classes ''CSF ; GM ; WM'' -bias_classes ''GM ; WM'' -bias_fwhm 25 -mixfrac UNI -main_N 5 -blur_meth BFT -prefix ' basePath];
    eval(str);
    
    % coregister segmentation with talairach
    origPath=pwd;cd(basePath);
    str='!@auto_tlrc -apar anat+tlrc -input Classes+orig.BRIK';
    eval(str);cd(origPath);
    
    
    % resample talairach segmentation on bold grid
    origPath=pwd;cd(basePath);
    str='!3dresample -master dsbold_e2_tlrc_al+tlrc -prefix resampled_Classes -input Classes+tlrc';
    eval(str);
    eval(str);cd(origPath);

end