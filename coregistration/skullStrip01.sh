#!/bin/bash

MPRAGE=$1
PD=$2
res=$3
coregFlag=$4

if [ -z "$1" ]
 then
  echo 'skull strip using PD. Inputs:'
  echo 'MPRAGE=$1, T1 filename' 
  echo 'PD=$2, PD filename'
  echo 'res=$3, output resolution (mm)'
  echo 'coregFlag=$4, coregister PD to MPRAGE? 0, 1'
  exit 1
fi


3dresample -master ${MPRAGE} -dxyz $res $res $res -inset ${MPRAGE} -prefix MPRAGE+orig

3dresample -master MPRAGE+orig -inset ${PD} -prefix PD+orig
	
	
if (( $coregFlag==1 ))
	then
	
	3dUnifize -prefix MPRAGE_uniform+orig MPRAGE+orig

	3dUnifize -prefix PD_uniform+orig PD+orig

	align_epi_anat.py -dset1 PD_uniform+orig \
			  -dset2 MPRAGE_uniform+orig \
			  -dset1to2 \
			  -dset1_strip 3dSkullStrip \
			  -dset2_strip 3dSkullStrip \
			  -cost lpc \
			  -Allineate_opts \
			  -final linear \
			  -warp shift_rotate \
			  -suffix _al \
			  -child_dset1 PD+orig
fi

if (( $coregFlag==0 ))
	then

	3dUnifize -prefix PD_uniform+orig PD+orig
	#3dcopy PD+orig PD_uniform+orig

	3dSkullStrip -input PD_uniform+orig -prefix PD_uniform_al+orig

fi


3dcalc -a MPRAGE+orig -b PD_uniform_al+orig -expr 'a*ispositive(b)' -prefix MPRAGE_ss+orig

3dAFNItoNIFTI MPRAGE_ss+orig

rm *.BRIK

rm *.HEAD

rm *.1D

