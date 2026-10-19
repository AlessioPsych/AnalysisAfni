to do:
arrange and clean up the code to run from a separate folder;
add the mrtrix processing (first step);
run the analysis on a surface based (along the visual hierarchy per ROI), or based on MNI space

to run this code you need to use a terminal open of the 'Flanker' folder,
a folder containing all the individual sub-xx data and the 'derivatives' folder.


01 run code forLoop.sh with the following command:
	tcsh forLoop.sh

this code runs the main processing across all participants (sub-01 to sub-26), calling the script 'sub_xx_afni_proc.sh' for each participant.
results are stored in the folder 'derivatives/processing_afni/'

02 run code forLoop_checkProcessing.sh with the following command:
	tcsh forLoop_checkProcessing.sh

this code runs creates some figure for the anatomy and EPI data for each participant (sub-01 to sub-26), calling the script 'visualizeProcessed.sh' for each participant.
results are stored in the folder 'derivatives/images_processing_afni/'. These figures are useful to quickly check whether some basic processing in terms of volreg, aligh and normalization went wrong.

03 run code stat_3dTtest_command.sh with the following command:
	tcsh stat_3dTtest_command.sh
this code runs some stats on the incongruent - congruent comparison across participants using 3dtstat++ and 3dMEMA

04 runFreesurferSuma.R
	runs freesurfer and suma on the anatomy of each participant
	
05 run code forLoop_checkSegmentations.sh with the following command:
	tcsh forLoop_checkSegmentations.sh

this code runs creates some figure for the anatomy segmentation for each participant (sub-01 to sub-26), calling the script 'visualizeSegmentations.sh' for each participant.
results are stored in the folder 'derivatives/images_processing_images_processing_Freesurfer/'. These figures are useful to quickly check whether some basic segmentation gave problems.

