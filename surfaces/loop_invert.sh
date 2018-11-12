#!/bin/bash


TRANSFORMATION="Transforms/amplitudeAnat_slicer_reg_al_reg_mat.aff12.1D"
ORIGINAL_SPACE="Aligned_to_anat/Dijk_cropped_anatomy.nii.gz"
FLAG_COMP_INV="0"
INTERP="NN"

if (($FLAG_COMP_INV=="1"))
then
	cat_matvec -ONELINE $TRANSFORMATION -I > tmp_inv.1D
else
	cp $TRANSFORMATION tmp_inv.1D	 
fi

for var in "$@"
do
	TRANSFORMED_VOL=$var
	echo $TRANSFORMED_VOL
	OUTPUT_VOL=${TRANSFORMED_VOL%%.*}
	OUTPUT_VOL=${OUTPUT_VOL%%+*}
	OUTPUT_VOL+="_inv"
	echo $OUTPUT_VOL
	
	if [ ! -f $OUTPUT_VOL".nii.gz" ]; then
   
    echo "File not found,creating"
    
	3dAllineate -1Dmatrix_apply "tmp_inv.1D" -master $ORIGINAL_SPACE -input $TRANSFORMED_VOL -prefix $OUTPUT_VOL".nii.gz" -interp $INTERP -final $INTERP
	
	else
	echo "File found, skipping"
	fi
done



