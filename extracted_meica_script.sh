#!/bin/sh

#  extracted_meica_script.sh
#  
#
#  Created by jacek on 11/17/17.
#

##########################################
3dWarp -overwrite -prefix /Users/jacek/Documents/MATLAB/LLLT/data/S06/NII/mebold2go/anat_do.nii.gz -deoblique /Users/jacek/Documents/MATLAB/LLLT/data/S06/NII/mebold2go/anat.nii;

# anat_u is in alignment with anat_do (?)
3dUnifize -overwrite -prefix /Users/jacek/Documents/MATLAB/LLLT/data/S06/NII/mebold2go/anat_u.nii.gz /Users/jacek/Documents/MATLAB/LLLT/data/S06/NII/mebold2go/anat_do.nii.gz;

# anat_ns is in alignment with anat_do
3dSkullStrip  -shrink_fac_bot_lim 0.3 -orig_vol -overwrite -prefix /Users/jacek/Documents/MATLAB/LLLT/data/S06/NII/mebold2go/anat_ns.nii.gz -input /Users/jacek/Documents/MATLAB/LLLT/data/S06/NII/mebold2go/anat_u.nii.gz;


# image is cropped
3dAutobox -overwrite -prefix /Users/jacek/Documents/MATLAB/LLLT/data/S06/NII/mebold2go/anat_ns.nii.gz /Users/jacek/Documents/MATLAB/LLLT/data/S06/NII/mebold2go/anat_ns.nii.gz;

## note: final anatomical is anat_ns.nii.gz


##########################################
# functional processing begins
# this just copies the bold files over to /meica.bold_e123/
3dcalc -a /Users/jacek/Documents/MATLAB/LLLT/data/S06/NII/mebold2go/bold_e1.nii -expr 'a' -prefix ./bold_e1.nii

nifti_tool -mod_hdr -mod_field sform_code 1 -mod_field qform_code 1 -infiles ./bold_e1.nii -overwrite

3dcalc -a /Users/jacek/Documents/MATLAB/LLLT/data/S06/NII/mebold2go/bold_e2.nii -expr 'a' -prefix ./bold_e2.nii

nifti_tool -mod_hdr -mod_field sform_code 1 -mod_field qform_code 1 -infiles ./bold_e2.nii -overwrite

3dcalc -a /Users/jacek/Documents/MATLAB/LLLT/data/S06/NII/mebold2go/bold_e3.nii -expr 'a' -prefix ./bold_e3.nii

nifti_tool -mod_hdr -mod_field sform_code 1 -mod_field qform_code 1 -infiles ./bold_e3.nii -overwrite

# take the latest anatomical and make it oblique like the bold
3dWarp -verb -card2oblique ./bold_e1.nii[0] -overwrite  -newgrid 1.000000 -prefix ./anat_ob.nii.gz /Users/jacek/Documents/MATLAB/LLLT/data/S06/NII/mebold2go/anat_ns.nii.gz

# copy some transformation matrix to disk
\grep  -A 4 '# mat44 Obliquity Transformation ::'  > bold_e_obla2e_mat.1D

##########################################
## this block for calculating and saving motion parameters
# remove impulses from bold_e1 (move along)
3dDespike -overwrite -prefix ./bold_e1_vrA.nii.gz ./bold_e1.nii

# data bricks oriented as axial slices (outdated function)
3daxialize -overwrite -prefix ./bold_e1_vrA.nii.gz ./bold_e1_vrA.nii.gz

# copy first frame of deimpulsed and axialized bold_e1 to disk
3dcalc -a ./bold_e1_vrA.nii.gz[0]  -expr 'a' -prefix eBbase.nii.gz

# align all bold frames to the first bold frame, which is saved in eBbase.nii.gz; save the results to a file
3dvolreg -overwrite -tshift -quintic  -prefix ./bold_e1_vrA.nii.gz -base eBbase.nii.gz -dfile ./bold_e1_vrA.1D -1Dmatrix_save ./bold_e_vrmat.aff12.1D ./bold_e1_vrA.nii.gz
1dcat './bold_e1_vrA.1D[1..6]{0..$}' > motion.1D

##########################################
## this block begins functional processing of echo 1
# despike bold_e1 and save it as bold_e1_pt
3dDespike -overwrite -prefix ./bold_e1_pt.nii.gz bold_e1.nii

# slice time alignment and store result as e1_ts+orig
3dTshift -heptic  -prefix ./e1_ts+orig ./bold_e1_pt.nii.gz

# fixes errors in the header
3drefit -view orig e1_ts*HEAD

# makes the slices axial
3daxialize  -overwrite -prefix ./e1_ts+orig ./e1_ts+orig

# removes any obliqueness in the e1_ts
3drefit -deoblique -TR 2.8 e1_ts+orig

##########################################
# repeats the above block for echo 2
3dDespike -overwrite -prefix ./bold_e2_pt.nii.gz bold_e2.nii
3dTshift -heptic  -prefix ./e2_ts+orig ./bold_e2_pt.nii.gz
3drefit -view orig e2_ts*HEAD
3daxialize  -overwrite -prefix ./e2_ts+orig ./e2_ts+orig
3drefit -deoblique -TR 2.8 e2_ts+orig

##########################################
# repeats the above block for echo 3
3dDespike -overwrite -prefix ./bold_e3_pt.nii.gz bold_e3.nii
3dTshift -heptic  -prefix ./e3_ts+orig ./bold_e3_pt.nii.gz
3drefit -view orig e3_ts*HEAD
3daxialize  -overwrite -prefix ./e3_ts+orig ./e3_ts+orig
3drefit -deoblique -TR 2.8 e3_ts+orig

##########################################
# prepares T2* and S0 volumes for use in functional masking and anat-func coregistration

# motion correction of the preprocessed bolds
# takes the functional in e1_ts+orig and aligns it with the base eBbase
3dAllineate -overwrite -final NN -NN -float -1Dmatrix_apply bold_e_vrmat.aff12.1D'{0..5}' -base eBbase.nii.gz -input e1_ts+orig'[0..5]' -prefix e1_vrA.nii.gz
3dAllineate -overwrite -final NN -NN -float -1Dmatrix_apply bold_e_vrmat.aff12.1D'{0..5}' -base eBbase.nii.gz -input e2_ts+orig'[0..5]' -prefix e2_vrA.nii.gz
3dAllineate -overwrite -final NN -NN -float -1Dmatrix_apply bold_e_vrmat.aff12.1D'{0..5}' -base eBbase.nii.gz -input e3_ts+orig'[0..5]' -prefix e3_vrA.nii.gz

# concatenates the 3 motion-corrected bold volumes along the Z-dimension and saves as basestack.nii.gz
3dZcat -prefix basestack.nii.gz  e1_vrA.nii.gz e2_vrA.nii.gz e3_vrA.nii.gz

# compute the T2s map for the concatenated bolds; must produce something called ocv.nii
/usr/bin/python /Users/jacek/abin/meica.libs/t2smap.py -d basestack.nii.gz -e 13,34,55

# do some operations on ocv.nii -- optimally combined volume
3dUnifize -prefix ./ocv_uni+orig ocv.nii

# skull strip optimally combined volume
3dSkullStrip -prefix ./ocv_ss.nii.gz -overwrite -input ocv_uni+orig

# I think this applies the skull stripping to t2svm and s0v
3dcalc -overwrite -a t2svm.nii -b ocv_ss.nii.gz -expr 'a*ispositive(a)*step(b)' -prefix t2svm_ss.nii.gz
3dcalc -overwrite -a s0v.nii -b ocv_ss.nii.gz -expr 'a*ispositive(a)*step(b)' -prefix s0v_ss.nii.gz
3daxialize -overwrite -prefix t2svm_ss.nii.gz t2svm_ss.nii.gz
3daxialize -overwrite -prefix ocv_ss.nii.gz ocv_ss.nii.gz
3daxialize -overwrite -prefix s0v_ss.nii.gz s0v_ss.nii.gz


##########################################
# align anat_ns to an atlas and create anat_ns_at
\@auto_tlrc -no_ss -init_xform AUTO_CENTER -base ${templateloc}/MNI_caez_N27+tlrc -input anat_ns.nii.gz -suffix _at

#
3dcopy /Users/jacek/Documents/MATLAB/LLLT/data/S06/NII/mebold2go/anat_ns_at.nii.gz anat_ns_at

# changes the "view" to be "orig" ?
3drefit -view orig anat_ns_at+tlrc
3dAutobox -prefix ./abtemplate.nii.gz ${templateloc}/MNI_caez_N27+tlrc
echo --------"Using alignp_mepi_anat.py to drive T2*-map weighted anatomical-functional coregistration"

# axializes the anatomical with the oblique matching the bold
3daxialize -overwrite -prefix ./anat_ob.nii.gz ./anat_ob.nii.gz

/usr/bin/python /Users/jacek/abin/meica.libs/alignp_mepi_anat.py -t t2svm_ss.nii.gz -a anat_ob.nii.gz -p mepi

cp alignp.mepi/mepi_al_mat.aff12.1D ./anat_al_mat.aff12.1D

# just the anatomical ?
cat_matvec -ONELINE /Users/jacek/Documents/MATLAB/LLLT/data/S06/NII/mebold2go/anat_ns_at.nii.gz::WARP_DATA -I > /Users/jacek/Documents/MATLAB/LLLT/data/S06/NII/mebold2go/anat_ns2at.aff12.1D

#
cat_matvec -ONELINE  /Users/jacek/Documents/MATLAB/LLLT/data/S06/NII/mebold2go/anat_ns_at.nii.gz::WARP_DATA -I bold_e_obla2e_mat.1D anat_al_mat.aff12.1D -I > bold_e_wmat.aff12.1D

cat_matvec -ONELINE  /Users/jacek/Documents/MATLAB/LLLT/data/S06/NII/mebold2go/anat_ns_at.nii.gz::WARP_DATA -I bold_e_obla2e_mat.1D anat_al_mat.aff12.1D -I  bold_e_vrmat.aff12.1D  > bold_e_vrwmat.aff12.1D





