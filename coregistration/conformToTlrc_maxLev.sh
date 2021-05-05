#!/bin/bash
# works with skull stripped data,
# see preProcessingPd.sh for that

MPRAGE=$1
DIR=$2
FLAGNONLIN=$3
MAXLEV=$4

if [ -z "$1" ]
 then
  echo 'Bash script to coregister a T1 volume on the TT_N27+tlrc. Inputs:'
  echo 'MPRAGE=$1, T1 filename' 
  echo 'DIR=$2, path_to/TT_N27+tlrc or path_to/TT_icbm452+tlrc ... TT_avg152T1+tlrc TT_EPI+tlrc'
  echo 'FLAGNONLIN=$3, perform non linear alignment'
  echo 'MAXLEV=$4, maxlev argument'
  exit 1
fi

3dcopy $MPRAGE MPRAGE.nii

atlasPath=$(printf '%s' $DIR) # set atlas file
matName=$(printf '%s_at.Xat.1D' ${MPRAGE%.*})
outname=$(printf '%s_at.nii' ${MPRAGE%.*})

echo $matName
echo $atlasPath
#echo $outname

@auto_tlrc -base $atlasPath -input MPRAGE.nii -no_ss -init_xform AUTO_CENTER

#mv $matName MPRAGE_tlrc_transform.Xat.1D
rm MPRAGE.nii

if (( $FLAGNONLIN==1 ))
 then
#	3dQwarp -blur 0 0 -iwarp -maxlev 3 \
#            -base $atlasPath -source MPRAGE_at.nii
	3dQwarp -iwarp -blur 0 0 -maxlev $MAXLEV \
            -base $atlasPath -source MPRAGE_at.nii

fi




