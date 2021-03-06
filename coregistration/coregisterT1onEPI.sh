#!/bin/bash

MPRAGE=$1 #'occipital_pole.nii.gz'
EPI=$2 #'amplitudeAnatomy.nii'
SEG=$3 #'occipital_pole_seg.nii.gz'
LEFTWM=$4
RIGHTWM=$5
RUNBBR=$6
CLIP=$7
CLIPDIRECTION=$8

3dcopy ${MPRAGE} MPRAGE+orig
3dcopy ${EPI} EPI+orig
3dcopy ${SEG} SEG+orig

instr=$(printf '@clip_volume -%s %d -input MPRAGE+orig -verb -prefix MPRAGE_crop+orig' $CLIPDIRECTION $CLIP )
$instr
instr=$(printf '@clip_volume -%s %d -input SEG+orig -verb -prefix SEG_crop+orig' $CLIPDIRECTION $CLIP )
$instr

3dZeropad -z 10 -prefix MPRAGE_zp+orig MPRAGE_crop+orig
3dZeropad -z 10 -prefix EPI_zp+orig EPI+orig
3dZeropad -z 10 -prefix SEG_zp+orig SEG_crop+orig

@Align_Centers -base EPI_zp+orig -cm -dset MPRAGE_zp+orig
cp MPRAGE_zp_shft+orig.HEAD SEG_zp_shft+orig.HEAD
cp SEG_zp+orig.BRIK SEG_zp_shft+orig.BRIK

3dcalc -a MPRAGE_zp_shft+orig -b SEG_zp_shft+orig -expr 'a*step(b)' -prefix MPRAGE_zp_shft_mask+orig

instr=$(printf '3dcalc -a SEG_zp_shft+orig -expr \u0027or( within(a,%d,%d), within(a,%d,%d) )\u0027 -prefix WM_zp_shft_mask+orig' $LEFTWM $LEFTWM $RIGHTWM $RIGHTWM )
Rscript $AFNI_TOOLBOXDIR/runCommandViaR.R $instr

3dAutomask -apply_prefix EPI_zp_ss+orig -clfrac 0.5 EPI_zp+orig

3dUniformize -anat MPRAGE_zp_shft_mask+orig -prefix MPRAGE_zp_shft_mask_unif+orig
3dUniformize -anat EPI_zp_ss+orig -prefix EPI_zp_ss_unif+orig


align_epi_anat.py -anat MPRAGE_zp_shft_mask_unif+orig -epi EPI_zp_ss_unif+orig -epi_base 0 -anat_has_skull no -epi_strip None -big_move -cost lpc

#align_epi_anat.py -anat MPRAGE_zp_shft_mask_unif+orig -epi EPI_zp_ss_unif+orig -epi_base 0 -anat_has_skull no -epi_strip None -big_move -volreg_opts -final wsinc5 -cost lpc

3dAllineate -final NN -1Dmatrix_apply MPRAGE_zp_shft_mask_unif_al_mat.aff12.1D -prefix WM_zp_shft_mask_al+orig WM_zp_shft_mask+orig

3dAllineate -final NN -1Dmatrix_apply MPRAGE_zp_shft_mask_unif_al_mat.aff12.1D -prefix SEG_zp_shft_al+orig SEG_zp_shft+orig




if (( $RUNBBR == 1 )) 

	then

	3dresample -master MPRAGE_zp_shft_mask_unif_al+orig -prefix EPI_zp_ss_unif_upsample+orig -inset EPI_zp_ss_unif+orig
	3dAFNItoNIFTI MPRAGE_zp_shft_mask_unif_al+orig
	3dAFNItoNIFTI EPI_zp_ss_unif_upsample+orig
	3dAFNItoNIFTI WM_zp_shft_mask_al+orig
#	3dAFNItoNIFTI DEPTH_zp_shft_al+orig
	3dAFNItoNIFTI SEG_zp_shft_al+orig

	echo "Running BBR"
	$FSLDIR/bin/flirt -ref MPRAGE_zp_shft_mask_unif_al.nii -in EPI_zp_ss_unif_upsample.nii -dof 6 -cost bbr -wmseg WM_zp_shft_mask_al.nii -omat MPRAGE_tMat.1D -out registeredEpi.nii -schedule ${FSLDIR}/etc/flirtsch/bbr.sch

	gunzip registeredEpi.nii

	convert_xfm -omat MPRAGE_tMat_inv.1D -inverse MPRAGE_tMat.1D

	$FSLDIR/bin/flirt -in MPRAGE_zp_shft_mask_unif_al.nii -ref MPRAGE_zp_shft_mask_unif_al.nii -applyxfm -init MPRAGE_tMat_inv.1D -interp nearestneighbour -out MPRAGE_zp_shft_mask_unif_al_bbr.nii

#	$FSLDIR/bin/flirt -in DEPTH_zp_shft_al.nii -ref MPRAGE_zp_shft_mask_unif_al.nii -applyxfm -init tMat_inv.1D -interp nearestneighbour -out DEPTH_zp_shft_al_bbr.nii

	$FSLDIR/bin/flirt -in SEG_zp_shft_al.nii -ref MPRAGE_zp_shft_mask_unif_al.nii -applyxfm -init tMat_inv.1D -interp nearestneighbour -out SEG_zp_shft_al_bbr.nii

#	$FSLDIR/bin/applywarp -i MPRAGE_zp_shft_mask_unif_al.nii -r MPRAGE_zp_shft_mask_unif_al.nii -o MPRAGE_zp_shft_mask_unif_al_bbr.nii --interp=nn --premat=tMat.1D

#	$FSLDIR/bin/applywarp -i DEPTH_zp_shft_al.nii -r MPRAGE_zp_shft_mask_unif_al.nii -o DEPTH_zp_shft_al_bbr.nii --interp=nn --premat=tMat.1D

#	$FSLDIR/bin/applywarp -i SEG_zp_shft_al.nii -r MPRAGE_zp_shft_mask_unif_al.nii -o SEG_zp_shft_al_bbr.nii --interp=nn --premat=tMat.1D

	echo 'interpolating to the EPI'
	3dresample -inset MPRAGE_zp_shft_mask_unif_al_bbr.nii.gz -prefix MPRAGE_zp_shft_mask_unif_al_bbr+orig -master EPI_zp_ss_unif+orig -rmode NN
#	3dresample -inset DEPTH_zp_shft_al_bbr.nii.gz -prefix DEPTH_zp_shft_al_bbr+orig -master EPI_zp_ss_unif+orig -rmode NN
	3dresample -inset SEG_zp_shft_al_bbr.nii.gz -prefix SEG_zp_shft_al_bbr+orig -master EPI_zp_ss_unif+orig -rmode NN

	3dZeropad -z -10 -prefix MPRAGE_zp_shft_mask_unif_al_bbr_epi+orig MPRAGE_zp_shft_mask_unif_al_bbr+orig
#	3dZeropad -z -10 -prefix DEPTH_zp_shft_al_bbr_epi+orig DEPTH_zp_shft_al_bbr+orig
	3dZeropad -z -10 -prefix SEG_zp_shft_al_bbr_epi+orig SEG_zp_shft_al_bbr+orig

	3dAFNItoNIFTI -prefix MPRAGE_zp_shft_mask_unif_al_bbr_epi MPRAGE_zp_shft_mask_unif_al_bbr_epi+orig 
#	3dAFNItoNIFTI -prefix DEPTH_zp_shft_al_bbr_epi DEPTH_zp_shft_al_bbr_epi+orig 
	3dAFNItoNIFTI -prefix SEG_zp_shft_al_bbr_epi SEG_zp_shft_al_bbr_epi+orig

	gzip MPRAGE_zp_shft_mask_unif_al_bbr_epi.nii
	gzip SEG_zp_shft_al_bbr_epi.nii

	rm MPRAGE_zp_shft_mask_unif_al_bbr.nii.gz
	rm DEPTH_zp_shft_al_bbr.nii.gz
	rm SEG_zp_shft_al_bbr.nii.gz
	rm registeredEpi.nii
	rm MPRAGE_zp_shft_mask_unif_al.nii
	rm EPI_zp_ss_unif_upsample.nii
	rm WM_zp_shft_mask_al.nii

	rm SEG_zp_shft_al.nii


	rm p*.1D
	rm v*.1D

	

else

	echo 'interpolating to the EPI'
	3dresample -inset MPRAGE_zp_shft_mask_unif_al+orig -prefix MPRAGE_zp_shft_mask_unif_al_resample+orig -master EPI_zp_ss_unif+orig -rmode NN

	3dresample -inset SEG_zp_shft_al+orig -prefix SEG_zp_shft_al_resample+orig -master EPI_zp_ss_unif+orig -rmode NN

	3dZeropad -z -10 -prefix MPRAGE_zp_shft_mask_unif_al_epi+orig MPRAGE_zp_shft_mask_unif_al_resample+orig

	3dZeropad -z -10 -prefix SEG_zp_shft_al_epi+orig SEG_zp_shft_al_resample+orig

	3dAFNItoNIFTI -prefix MPRAGE_al_epi MPRAGE_zp_shft_mask_unif_al_epi+orig 

	3dAFNItoNIFTI -prefix SEG_al_epi SEG_zp_shft_al_epi+orig

	gzip MPRAGE_al_epi.nii

	gzip SEG_al_epi.nii

	rm MPRAGE*.BRIK
	rm MPRAGE*.HEAD
	rm EPI*.BRIK
	rm EPI*.HEAD
	rm SEG*.BRIK
	rm SEG*.HEAD
	rm WM*.BRIK
	rm WM*.HEAD
	rm p*.1D
	rm v*.1D


fi

