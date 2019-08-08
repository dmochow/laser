% using this to clean up data structure
clear all; close all; clc
%subjStrs={'S02','S03','S04','S05','S06','S07','S09','S10','S11','S12','S13','S14','S15','S16','S17','S18','S19','S20','S21','S22'};
subjStrs={'S23','S24','S25','S26','S27','S28','S29','S30'};
sourceDir='/Users/jacekdmochowski/PROJECTS/LLLT/data/';
targetDir='/Users/jacekdmochowski/PROJECTS/LLLT/data/scanner_data/';
nSubjects=numel(subjStrs);
for s=1:nSubjects
    subjStr=subjStrs{s};
    mkdir([targetDir subjStr]);
    
    
%     copyfile([sourceDir subjStr '/NII/anat.nii'],[targetDir subjStr '/anat.nii']);
%     copyfile([sourceDir subjStr '/NII/bold_e1.nii'],[targetDir subjStr '/bold_e1.nii']);
%     copyfile([sourceDir subjStr '/NII/bold_e2.nii'],[targetDir subjStr '/bold_e2.nii']);
%     copyfile([sourceDir subjStr '/NII/bold_e3.nii'],[targetDir subjStr '/bold_e3.nii']);
%     copyfile([sourceDir subjStr '/BIOPAC/biopac.mat'],[targetDir subjStr '/biopac.mat']);
%     delete([sourceDir subjStr '/NII/*.*']);
    
    copyfile([sourceDir subjStr '/anat.nii'],[targetDir subjStr '/anat.nii']);
    copyfile([sourceDir subjStr '/bold_e1.nii'],[targetDir subjStr '/bold_e1.nii']);
    copyfile([sourceDir subjStr '/bold_e2.nii'],[targetDir subjStr '/bold_e2.nii']);
    copyfile([sourceDir subjStr '/bold_e3.nii'],[targetDir subjStr '/bold_e3.nii']);
    copyfile([sourceDir subjStr '/*.mat'],[targetDir subjStr '/']);
    delete([sourceDir subjStr '/*.*']);
end