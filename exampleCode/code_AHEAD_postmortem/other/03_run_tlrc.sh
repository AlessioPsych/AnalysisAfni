#!/bin/bash

# to run, remember to change the participant folder in line 12 to either 
#122017_resample/ 
#or
#152017_resample/

mainFolder="/media/alessiofracasso/DATADRIVE1/AHEAD_exvivo/"
atlasFolder="suma_MNI152_2009/"
atlasFile="brain.nii"
atlasRois="aparc+aseg.nii"
dataFolder="152017_resample"

cd $mainFolder
cd $dataFolder

if [ -f $atlasFile ]; then
	echo "removing file $atlasFile" 
	rm $atlasFile
fi
cp $mainFolder/$atlasFolder/$atlasFile $mainFolder/$dataFolder/$atlasFile 

if [ -f $atlasRois ]; then
	echo "removing file $atlasRois" 
	rm $atlasRois
fi
cp $mainFolder/$atlasFolder/$atlasRois $mainFolder/$dataFolder/$atlasRois 

#@auto_tlrc -init_xform CENTER_CM -no_ss -base $atlasFile -input qR1_resampled.nii.gz
@auto_tlrc -init_xform CENTER_CM -no_ss -base $atlasFile -input PD_resampled.nii.gz

#3dQwarp -iwarp -blur 0 0 -maxlev 1 \
#            -base $atlasFile -source PD_resampled_at.nii.gz
            
#3dNwarpApply 
            
3dcalc -a $atlasRois -expr 'within(a,1034.5,1035.5)' -prefix insulaLeft.nii.gz
3dcalc -a $atlasRois -expr 'within(a,2034.5,2035.5)' -prefix insulaRight.nii.gz

cat_matvec PD_resampled_at.Xat.1D -I > atlasToVolume_mat.1D 

3dAllineate -base PD_resampled.nii.gz -input insulaLeft.nii.gz -prefix insulaLeft_transformed.nii.gz -1Dmatrix_apply atlasToVolume_mat.1D -final NN

3dAllineate -base PD_resampled.nii.gz -input insulaRight.nii.gz -prefix insulaRight_transformed.nii.gz -1Dmatrix_apply atlasToVolume_mat.1D -final NN

