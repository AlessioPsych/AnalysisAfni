#!/bin/bash

# to run, run individually for each participant: 
#122017_resample/ from line 8 to 21
#or
#152017_resample/ from line 25 to 36

mainFolder="/media/alessiofracasso/DATADRIVE1/AHEAD_exvivo/"
dataFolder="152017_resample"

cd $mainFolder
cd $dataFolder

# cruise, gray matter = 1, white matter = 2 

3dcalc -a qR1_resampled.nii.gz -b cruiseLeft_resampled.nii.gz -c cruiseRight_resampled.nii.gz -expr 'a*0.9*within(b,0.5,1.5) + a*0.9*within(c,0.5,1.5) + a*1.2*within(b,1.5,2.5) + a*1.2*within(c,1.5,2.5) + a*not(step(b+c))' -prefix qR1_resampled_inputFreesurfer.nii.gz

SUBJECTS_DIR="/media/alessiofracasso/DATADRIVE1/AHEAD_exvivo/152017_resample/"
echo $SUBJECTS_DIR
recon-all -subjid Freesurfer_result -i qR1_resampled_inputFreesurfer.nii.gz -all -parallel -openmp 8
@SUMA_Make_Spec_FS -NIFTI -fspath $SUBJECTS_DIR/Freesurfer_result -sid Freesurfer_result



mainFolder="/media/alessiofracasso/DATADRIVE1/AHEAD_exvivo/"
dataFolder="122017_resample"

3dcalc -a qR1_resampled.nii.gz -b cruiseLeft_resampled.nii.gz -c cruiseRight_resampled.nii.gz -expr 'a*0.8*within(b,0.5,1.5) + a*0.8*within(c,0.5,1.5) + a*1.3*within(b,1.5,2.5) + a*1.3*within(c,1.5,2.5) + a*not(step(b+c))' -prefix qR1_resampled_inputFreesurfer.nii.gz

3dUnifize -prefix qR1_resampled_inputFreesurfer_uni.nii.gz qR1_resampled_inputFreesurfer.nii.gz

SUBJECTS_DIR="/media/alessiofracasso/DATADRIVE1/AHEAD_exvivo/122017_resample/"
echo $SUBJECTS_DIR
recon-all -subjid Freesurfer_result -i qR1_resampled_inputFreesurfer_uni.nii.gz -bigventricles -all -parallel -openmp 8
# -nuiterations 6
@SUMA_Make_Spec_FS -NIFTI -fspath $SUBJECTS_DIR/Freesurfer_result -sid Freesurfer_result


