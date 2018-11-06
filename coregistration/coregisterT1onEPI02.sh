#!/bin/bash

MPRAGE=$1 
EPI=$2 
ALIGNCENTERS=$3
COSTFUN=$4
FLAGNONLINEAR=$5
UNIF=$6
PATCHMIN=$7
FLAGREMOVEINTERMEDIATE=$8

if [ -z "$1" ]
 then
  echo 'Bash script to coregister a T1 volume on an EPI using affine + non linear transforms. Inputs:'
  echo 'MPRAGE=$1, T1 filename' 
  echo 'EPI=$2, EPI filename'
  echo 'ALIGNCENTERS=$3, [1 center of volume; 2 center of mass]'
  echo 'COSTFUN=$4, cost function: lpc, nmi, mi etc etc (it uses function 3dAllineate)'
  echo 'FLAGNONLINEAR=$5, perform non-linear coregistration [0=no, 1=yes, qwarp, 2= yes, qwarp plusminus]'
  echo 'UNIF=$6, make transformed T1 volume uniform [0 (no), 1 (yes, 3dUnifize), 2 (yes, 3dUnifize inverted contrast), 3 (yes, 3dUniformize) ]'
  echo 'PATCHMIN=$7, minimum patch size (mm) for 3dQwarp'
  echo 'FLAGREMOVEINTERMEDIATE=$8, remove temp datasets? [0, 1]'
  exit 1
fi

3dcopy $MPRAGE coreg_tempFile_MPRAGE+orig
3dcopy $EPI coreg_tempFile_EPI+orig

if (( $ALIGNCENTERS==1 ))
	then
	@Align_Centers -base coreg_tempFile_EPI+orig -dset coreg_tempFile_MPRAGE+orig
fi
if (( $ALIGNCENTERS==2 ))
	then
	@Align_Centers -base coreg_tempFile_EPI+orig -cm -dset coreg_tempFile_MPRAGE+orig
fi

align_epi_anat.py -anat coreg_tempFile_MPRAGE_shft+orig -epi coreg_tempFile_EPI+orig -epi_base 0 -anat_has_skull no -epi_strip None -big_move -cost $COSTFUN

if (( $FLAGNONLINEAR==1 ))
	then

	preProcess02.sh coreg_tempFile_MPRAGE_shft_al+orig 0 0 0 0 0 coreg_tempFile_MPRAGE_T2Like+orig $UNIF 0


	3dresample -inset coreg_tempFile_MPRAGE_T2Like+orig -master coreg_tempFile_EPI+orig -prefix 	coreg_tempFile_MPRAGE_T2Like_resample+orig

	3dQwarp -base coreg_tempFile_EPI+orig -source coreg_tempFile_MPRAGE_T2Like_resample+orig 	-blur 0 0 -iwarp -patchmin $PATCHMIN

	cat_matvec -ONELINE coreg_tempFile_MPRAGE_shft_al_mat.aff12.1D coreg_tempFile_MPRAGE_shft.1D > combined.1D
cat_matvec -ONELINE combined.1D -I > combined_inv.1D

fi

if (( $FLAGNONLINEAR==2 ))
  then

	preProcess02.sh coreg_tempFile_MPRAGE_shft_al+orig 0 0 0 0 0 coreg_tempFile_MPRAGE_T2Like+orig $UNIF 0


	3dresample -inset coreg_tempFile_MPRAGE_T2Like+orig -master coreg_tempFile_EPI+orig -prefix 	coreg_tempFile_MPRAGE_T2Like_resample+orig

	3dQwarp -base coreg_tempFile_EPI+orig -source coreg_tempFile_MPRAGE_T2Like_resample+orig -plusminus	-blur 0 0 -iwarp -patchmin $PATCHMIN

	cat_matvec -ONELINE coreg_tempFile_MPRAGE_shft_al_mat.aff12.1D coreg_tempFile_MPRAGE_shft.1D > combined.1D
cat_matvec -ONELINE combined.1D -I > combined_inv.1D

fi


if (( $FLAGREMOVEINTERMEDIATE==1 ))
    then
        rm coreg_temp*.BRIK
        rm coreg_temp*.HEAD
        rm coreg_temp*.1D
fi

# check if 3dNwarp apply works also in this case (manual resampling, otherwise, use the -resample flag)
# check if the coregistration works on the oblique dataset (Spinoza data)

 

#3dNwarpApply -source anatomyClip_process.nii.gz -nwarp 'Qwarp_WARP+orig combined.1D' -prefix anatomySingleShot.nii.gz -master amplitudeAnatomyClip_process.nii.gz -interp wsinc5


#cat_matvec -ONELINE combined.1D -I > combined_inv.1D



#3dNwarpApply -source anatomySingleShot.nii.gz -nwarp 'combined_inv.1D Qwarp_WARPINV+orig' -prefix anatomySingleShot_inv.nii.gz -master anatomyClip_process.nii.gz -interp wsinc5

