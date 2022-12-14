#!/bin/bash

MPRAGE=$1
INV2=$2
RES=$3

if [ -z "$1" ]
then
echo 'segment script using nighres (MP2RAGE)'
echo 'Inputs:'
echo 'MPRAGE, input anatomy'
echo 'INV2, second inversion volume'
echo 'RES, output resolution, mm'
echo 'example call: skullStrip_mp2rage.sh mp2rage.nii.gz inv2.nii.gz'
exit 1
fi

3dWarp -deoblique -prefix _del_MP2RAGE_deob.nii.gz $MPRAGE
3dWarp -deoblique -prefix _del_MP2RAGE_INV2_deob.nii.gz $INV2

3dresample -dxyz $RES $RES $RES -orient RAI -prefix _del_MP2RAGE_RAI.nii.gz -inset _del_MP2RAGE_deob.nii.gz -rmode Lin
3dresample -master _del_MP2RAGE_RAI.nii.gz -prefix _del_INV2_RAI.nii.gz -inset _del_MP2RAGE_INV2_deob.nii.gz -rmode Lin
3dSkullStrip -input _del_INV2_RAI.nii.gz -prefix _del_SKULL_STRIP.nii.gz -shrink_fac 0.4
3dcalc -a _del_SKULL_STRIP.nii.gz -b _del_MP2RAGE_RAI.nii.gz -expr 'b*step(a)' -prefix _del_MP2RAGE_RAI_SKULL_STRIP.nii.gz 
3dAutobox -input _del_MP2RAGE_RAI_SKULL_STRIP.nii.gz -noclust -prefix MP2RAGE_RAI_SKULL_STRIP_BOX.nii.gz -npad 5 

3dcopy MP2RAGE_RAI_SKULL_STRIP_BOX.nii.gz anatCopy.nii.gz

rm _del_*
rm MP2RAGE_RAI_SKULL_STRIP_BOX.nii.gz


