#!/bin/bash

VOLUME=$1
ZEROPAD=$2
UNIF=$3
MASK=$4
CLFRAC=$5
DILATE=$6
OUTNAME=$7

3dcopy $VOLUME tempFile_vol+orig

3dZeropad -z $ZEROPAD -prefix tempFile_vol_clip_zp+orig tempFile_vol+orig

if (( $UNIF==0 ))
	then
	3dcopy tempFile_vol_clip_zp+orig tempFile_vol_clip_zp_unif+orig
fi

if (( $UNIF==1 ))
	then
	3dUnifize -input tempFile_vol_clip_zp+orig -prefix tempFile_vol_clip_zp_unif+orig
fi

if (( $UNIF==2 ))
	then
	3dUnifize -input tempFile_vol_clip_zp+orig -prefix tempFile_vol_clip_zp_unif+orig -T2 -T2
fi



if (( $MASK==0 ))
	then
	3dcopy tempFile_vol_clip_zp_unif+orig $OUTNAME
fi

if (( $MASK==1 ))
	then
	3dAutomask -clfrac $CLFRAC -dilate $DILATE -apply_prefix $OUTNAME tempFile_vol_clip_zp_unif+orig
fi

rm tempFile_*.BRIK
rm tempFile_*.HEAD

