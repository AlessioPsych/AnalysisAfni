# makes use of  3dWarp -disp_obl_xform_only -deoblique ttt_original_EPI_ts.nii  > ttt_mat_obl_transform.1D, it needs afni version above 2017 or so
# set directory
cd /analyse/Project0226/KastnerModel/staging_area/JFE27_08092022/coregister

# clean up
rm coregMat.1D
rm invMat.1D
rm singleShot.nii.gz
rm ttt*

cp /analyse/Project0226/KastnerModel/staging_area/JFE27_08092022/amplitudeAnatomy.nii amplitudeAnatomy.nii

mv amplitudeAnatomy.nii ttt_amplitudeAnatomy.nii

cp /analyse/Project0226/KastnerModel/staging_area/JFE27_08092022/EPI/07_GN19NE455_JFE27_08092022_MB_GRE_EPI_155VOL_BARS_MAPPING_50_slices_20220908140040_7.nii /analyse/Project0226/KastnerModel/staging_area/JFE27_08092022/coregister/ttt_original_EPI_ts.nii

# gets the oblique matrix and saves it on file as a 3X4 matrix
3dWarp -disp_obl_xform_only -deoblique ttt_original_EPI_ts.nii  > ttt_mat_obl_transform.1D

# removes the comments from the oblique matrix file
grep -o '^[^#]*' ttt_mat_obl_transform.1D > ttt_mat_obl_transform_no_comments.1D

# saves the oblique matrix file as a 1D file (as 3dAllineate requires)
cat_matvec -ONELINE ttt_mat_obl_transform_no_comments.1D -I > ttt_mat_obl_transform_no_comments_ONELINE.1D

# compute average original EPI, just to provide a space to test the deoblique matrix
3dTstat -mean -prefix ttt_original_EPI_mean.nii ttt_original_EPI_ts.nii

# deoblique average original EPI, just to provide a space to test the deoblique matrix
3dWarp -deoblique -prefix ttt_original_EPI_mean_deob.nii ttt_original_EPI_mean.nii

# apply transformation matrix to deoblique matrix, to check, here I am using the deoblique(d) space as master
3dAllineate -prefix ttt_original_EPI_mean_applyDeobMat.nii -1Dmatrix_apply ttt_mat_obl_transform_no_comments_ONELINE.1D -final wsinc5 -input ttt_original_EPI_mean.nii -master ttt_original_EPI_mean_deob.nii

# align centers of EPI and anatCopy.nii.gz
@Align_Centers -base anatCopy.nii.gz -dset ttt_amplitudeAnatomy.nii

# upsample EPI to resolution of ANATOMY
3dresample -prefix _ttt_epi_startingPoint_resample.nii -input ttt_amplitudeAnatomy_shft.nii -rmode Linear -master anatCopy.nii.gz

# apply deob matrix
3dAllineate -prefix _ttt_epi_startingPoint_resample_applyDeobMat.nii -1Dmatrix_apply ttt_mat_obl_transform_no_comments_ONELINE.1D -final wsinc5 -input _ttt_epi_startingPoint_resample.nii -master anatCopy.nii.gz


# ---!!! UPDATE THIS LINE !!!---
# apply nudge manually:
#	- afni, underlay, ttt_amplitudeAnatomy.nii
#	- new, _ttt_epi_startingPoint_resample.nii
# 	- define datamode, plugins, nudge dataset, choose dataset, _ttt_epi_startingPoint_resample
#	- print command line and copy-paste here
3drotate -quintic -clipit -rotate 0.00I 0.00R 0.00A -ashift -5.00S -6.00L -4.00P -prefix nudgedDataset.nii.gz _ttt_epi_startingPoint_resample_applyDeobMat.nii

# store affine transformation in 1D format
cat_matvec -ONELINE 'nudgedDataset.nii.gz::ROTATE_MATVEC_000000' -I > nudgeMat.1D
#cat_matvec -ONELINE 'amplitudeAnatomy_ts.nii::IJK_TO_DICOM_REAL' -I > obliqueMat.1D

#3dAllineate -prefix _ttt_epi_startingPoint_mask_shft_nudge_rotmat_inv.nii -1Dmatrix_apply nudgeMat.1D -final wsinc5 -input _ttt_epi_startingPoint_resample.nii -master _ttt_epi_startingPoint_resample.nii

# apply transformation matrix to starting point matrix
3dAllineate -prefix _ttt_epi_startingPoint_mask_shft_nudge_rotmat_inv.nii -1Dmatrix_apply nudgeMat.1D -final wsinc5 -input _ttt_epi_startingPoint_resample_applyDeobMat.nii -master _ttt_epi_startingPoint_resample_applyDeobMat.nii

# perform coregistration
align_epi_anat.py -epi_base 0 -epi _ttt_epi_startingPoint_mask_shft_nudge_rotmat_inv.nii -anat anatCopy.nii.gz -epi2anat -cost lpc -anat_has_skull no -epi_strip None -Allineate_opts -interp cubic -onepass -weight_frac 1.0 -maxrot 5 -maxshf 5 -VERB -warp aff

# store all transformations in single transform matrix
cat_matvec -ONELINE _ttt_epi_startingPoint_mask_shft_nudge_rotmat_inv_al_mat.aff12.1D nudgeMat.1D ttt_mat_obl_transform_no_comments_ONELINE.1D ttt_amplitudeAnatomy_shft.1D > coregMat.1D
cat_matvec -ONELINE coregMat.1D -I > invMat.1D

# apply single transform matrix to initial input
3dAllineate -prefix singleShot.nii.gz -1Dmatrix_apply coregMat.1D -final linear -input ttt_amplitudeAnatomy.nii -master anatCopy.nii.gz

# clean up
rm ttt*
rm _ttt*
rm nudgedDataset.nii.gz
rm nudgeMat.1D
rm anatCopy_al_mat.aff12.1D

