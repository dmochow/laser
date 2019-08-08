# bare bones pipeline
# 11/20/17
#
# Prantik's script is too much to figure out, so I am making my own script here and will try to build it up incrementally
#


#
align_epi_anat.py -anat /Users/jacek/Documents/MATLAB/LLLT/data/SDB_/anat.nii -epi /Users/jacek/Documents/MATLAB/LLLT/data/SDB_/bold_e2.nii -epi_base 5 -save_Al_in -save_vr -save_epi_ns

3dAllineate -cubic -1Dmatrix_apply /Users/jacek/Documents/MATLAB/LLLT/data/SDB_/anat_al_mat.aff12.1D -prefix /Users/jacek/Documents/MATLAB/LLLT/data/SDB_/roiMask_al /Users/jacek/Documents/MATLAB/LLLT/data/SDB_/roiMask.nii

