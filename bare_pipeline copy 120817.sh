# bare bones pipeline
# 11/22/17
#
# Prantik's script is too much to figure out, so I am making my own script here and will try to build it up incrementally
#

## first create tailarach registered anatomy

#@auto_tlrc -base TT_N27+tlrc -input anat.nii

#@auto_tlrc -apar anat_at.nii -input roiMask.nii

## now align the epis to the anatomy
#align_epi_anat.py -anat anat.nii -epi bold_e2.nii      \
#-epi_base 5 -child_epi bold_e1.nii bold_e3.nii  \
#-epi2anat -suffix _al_at  \
#-tlrc_apar anat_at.nii

## resample mask on the aligned and atlased bold
3dresample -master bold_e2_tlrc_al_at+tlrc -prefix roiMask_at_resampled -input roiMask_at.nii