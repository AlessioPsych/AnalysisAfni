#!/bin/bash

INPUTVOL=$1
CLFRAC=$2
DILATE=$3
ERODE=$4
ZEROPAD=$5
OUTVOLUME=$6

3dAutobox -input $INPUTVOL -prefix tempVol_autobox+orig

3dAutomask -clfrac $CLFRAC -dilate $DILATE -erode $ERODE -apply_prefix tempVol_autobox_masked+orig tempVol_autobox+orig

3dZeropad -R $ZEROPAD -L $ZEROPAD -A $ZEROPAD -P $ZEROPAD -I $ZEROPAD -S $ZEROPAD -prefix tempVol_autobox_masked_zp+orig tempVol_autobox_masked+orig

3dAFNItoNIFTI -prefix $OUTVOLUME tempVol_autobox_masked_zp+orig

rm tempVol*.BRIK
rm tempVol*.HEAD
