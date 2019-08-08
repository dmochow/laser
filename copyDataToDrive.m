clear all; close all; clc
subjStrs={'S04','S05','S07','S06','S09','S10','S11','S12','S13','S14','S15','S16','S17','S18','S19','S20','S21','S22'};
%subjStrs={'S06','S11'};
sourceDir='/Users/jacekdmochowski/PROJECTS/LLLT/data/';
targetDir='/Users/jacekdmochowski/Google Drive/Dmochowski Lab/Research Projects/Laser/data/';
nSubjects=numel(subjStrs);
for s=1:nSubjects
    subjStr=subjStrs{s};
    mkdir([targetDir subjStr]);
%     str=['!gunzip ' sourceDir subjStr  '/NII/anat+orig.BRIK.gz']; eval(str);
%     str=['!gunzip ' sourceDir subjStr  '/NII/brain_mask+tlrc.BRIK.gz']; eval(str);
%     str=['!gunzip ' sourceDir subjStr  '/NII/resampled_brain_mask+tlrc.BRIK.gz']; eval(str);
%     str=['!gunzip ' sourceDir subjStr  '/NII/resampled_roi_r19_z26+tlrc.BRIK.gz']; eval(str);
%     copyfile([sourceDir subjStr '/NII/anat+orig.HEAD'],[targetDir subjStr '/anat+orig.HEAD']);
%     copyfile([sourceDir subjStr '/NII/anat+orig.BRIK'],[targetDir subjStr '/anat+orig.BRIK']);
%     copyfile([sourceDir subjStr '/NII/anat+tlrc.HEAD'],[targetDir subjStr '/anat+tlrc.HEAD']);
%     copyfile([sourceDir subjStr '/NII/anat+tlrc.BRIK'],[targetDir subjStr '/anat+tlrc.BRIK']);
%     copyfile([sourceDir subjStr '/NII/brain_mask+orig.HEAD'],[targetDir subjStr '/brain_mask+orig.HEAD']);
%     copyfile([sourceDir subjStr '/NII/brain_mask+orig.BRIK'],[targetDir subjStr '/brain_mask+orig.BRIK']);
%     copyfile([sourceDir subjStr '/NII/brain_mask+tlrc.HEAD'],[targetDir subjStr '/brain_mask+tlrc.HEAD']);
%     copyfile([sourceDir subjStr '/NII/brain_mask+tlrc.BRIK'],[targetDir subjStr '/brain_mask+tlrc.BRIK']);
%     copyfile([sourceDir subjStr '/NII/resampled_brain_mask+tlrc.HEAD'],[targetDir subjStr '/resampled_brain_mask+tlrc.HEAD']);
%     copyfile([sourceDir subjStr '/NII/resampled_brain_mask+tlrc.BRIK'],[targetDir subjStr '/resampled_brain_mask+tlrc.BRIK']);
%     copyfile([sourceDir subjStr '/NII/roi_r19_z26+orig.HEAD'],[targetDir subjStr '/roi+orig.HEAD']);
%     copyfile([sourceDir subjStr '/NII/roi_r19_z26+orig.BRIK'],[targetDir subjStr '/roi+orig.BRIK']);
%     copyfile([sourceDir subjStr '/NII/roi_r19_z26+tlrc.HEAD'],[targetDir subjStr '/roi+tlrc.HEAD']);
%     copyfile([sourceDir subjStr '/NII/roi_r19_z26+tlrc.BRIK'],[targetDir subjStr '/roi+tlrc.BRIK']);
%     copyfile([sourceDir subjStr '/NII/resampled_roi_r19_z26+tlrc.HEAD'],[targetDir subjStr '/resampled_roi+tlrc.HEAD']);
%     copyfile([sourceDir subjStr '/NII/resampled_roi_r19_z26+tlrc.BRIK'],[targetDir subjStr '/resampled_roi+tlrc.BRIK']);
copyfile([sourceDir subjStr '/NII/oBoldEcho1+tlrc.HEAD'],[targetDir subjStr '/oBoldEcho1+tlrc.HEAD']);
copyfile([sourceDir subjStr '/NII/oBoldEcho1+tlrc.BRIK'],[targetDir subjStr '/oBoldEcho1+tlrc.BRIK']);
copyfile([sourceDir subjStr '/NII/oBoldEcho2+tlrc.HEAD'],[targetDir subjStr '/oBoldEcho2+tlrc.HEAD']);
copyfile([sourceDir subjStr '/NII/oBoldEcho2+tlrc.BRIK'],[targetDir subjStr '/oBoldEcho2+tlrc.BRIK']);
copyfile([sourceDir subjStr '/NII/oBoldEcho3+tlrc.HEAD'],[targetDir subjStr '/oBoldEcho3+tlrc.HEAD']);
copyfile([sourceDir subjStr '/NII/oBoldEcho3+tlrc.BRIK'],[targetDir subjStr '/oBoldEcho3+tlrc.BRIK']);
%     copyfile([sourceDir subjStr '/NII/resampled_anat+tlrc.BRIK'],[targetDir subjStr '/resampled_anat+tlrc.BRIK']);
%     copyfile([sourceDir subjStr '/NII/resampled_anat+tlrc.HEAD'],[targetDir subjStr '/resampled_anat+tlrc.HEAD']);
end