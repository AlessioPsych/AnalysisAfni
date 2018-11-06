#!/bin/bash

INPUTDIR=$1
OUTPUTDIR=$2

if [ -z "$1" ]
then
echo 'computes amplitude anatomy (mean of mean ts)'
echo 'reads .volreg files processed by motionCorrect.afni.sh'
echo 'Inputs:'
echo 'INPUTDIR=$1, input directory'
echo 'OUTPUTDIR=$2, output directory'
exit 1
fi


cd $INPUTDIR

mean_files=(**.volreg+orig**.BRIK)

for ((i=0; i<${#mean_files[@]}; i++)); do

   outputFileNameMean=$(printf '%s%d%s' 'meanAmplitude_' $i '+orig')
   echo ${outputFileNameMean}

   3dTstat -prefix ${outputFileNameMean} ${mean_files[i]}   

   3dAFNItoNIFTI ${outputFileNameMean}

   rm $(printf '%s%s' ${outputFileNameMean} '.HEAD')
   rm $(printf '%s%s' ${outputFileNameMean} '.BRIK')

   
done


mean_files_nii=( meanAmplitude_**.nii )

outputFileName=$(printf '3dMean -prefix amplitudeAnatomy+orig %s' ${mean_files_nii[0]} )

for ((i=1; i<=${#mean_files_nii[@]}; i++)); do
     
    outputFileName=$(printf '%s %s ' ${outputFileName} ${mean_files_nii[i]} )
   
done

echo ${outputFileName}

${outputFileName}

3dAFNItoNIFTI amplitudeAnatomy+orig

rm meanAmplitude**.nii
rm amplitudeAnatomy+orig.BRIK
rm amplitudeAnatomy+orig.HEAD

mv amplitudeAnatomy.nii ../$OUTPUTDIR

cd ..



