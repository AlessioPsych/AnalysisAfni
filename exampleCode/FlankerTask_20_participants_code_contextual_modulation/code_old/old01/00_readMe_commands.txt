to do:
arrange and clean up the code to run from a separate folder;
add the mrtrix processing (first step);
run the analysis on a surface based (along the visual hierarchy per ROI), or based on MNI space

01_00 run the code: 01_00_denoiseData.sh which denoises the data using dwidenoise from mrtrix, below is how to call the script from the mrtrix docker. The code copies the original functional data into derivatives/mrtrix3 and removes thermal noise
	# Instructions to run:
	sudo docker run -v /media/alessiofracasso/DATADRIVE1/Flanker:/mrtrixDataFolder --rm -it mrtrix3/mrtrix3
	cd /mrtrixDataFolder
	sh 01_00_denoiseData.sh

01_01 run the code: 01_01_copyAnatInDenoisedFolder.sh which copies the anatomy folder from the original data to the mrtrix3 denoised folder. In essence 01_00_denoiseData.sh and 01_01_copyAnatInDenoisedFolder.sh create a new denoised dataset, with the same folder structure as the original data.
	# Instructions to run, from the code/ folder:
	sh 01_01_copyAnatInDenoisedFolder.sh

02 run code 02_00_forLoop_MNI_denoised.sh with the command below. This code runs the main processing across all participants (sub-01 to sub-26). Results are stored in the folder 'derivatives/processing_afni_MNI_denoised/'. This code calls the code: 02_01_sub_xx_afni_proc_MNI_denoised.sh.
to run, from the main folder:
	# tcsh 02_00_forLoop_MNI_denoised.sh

03 run code 03_00_forLoop_MNI_no_denoised.sh with the command below. This code runs the main processing across all participants (sub-01 to sub-26). Results are stored in the folder 'derivatives/processing_afni_MNI_no_denoised/'. This code calls the code: 03_01_sub_xx_afni_proc_MNI_no_denoised.sh.
to run, from the main folder:
	# tcsh 03_00_forLoop_MNI_no_denoised.sh

04 run code 04_00_forLoop_denoised.sh with the command below. This code runs the main processing across all participants (sub-01 to sub-26). Results are stored in the folder 'derivatives/processing_afni_denoised/'. This code calls the code: 04_01_sub_xx_afni_proc_denoised.sh.
to run, from the main folder:
	# tcsh 04_01_forLoop_denoised.sh

05 run code 05_00_forLoop_no_denoised.sh with the command below. This code runs the main processing across all participants (sub-01 to sub-26). Results are stored in the folder 'derivatives/processing_afni_no_denoised/'. This code calls the code: 05_01_sub_xx_afni_proc_no_denoised.sh. To run, from the main folder:
	# tcsh 05_01_forLoop_no_denoised.sh
	
06 06_runFreesurferSuma.R runs freesurfer and suma on the anatomy of each participant
	Rscript 06_runFreesurferSuma.R
	
07 run code 07_00_forLoop_checkProcessing_MNI.sh this code runs creates some figure for the anatomy and EPI data for each participant (sub-01 to sub-26), calling the script '07_01_visualizeProcessed_MNI.sh' for each participant. Results are stored in the folder 'derivatives/images_processing_afni_MNI_denoised/'. These figures are useful to quickly check whether some basic processing in terms of volreg, aligh and normalization went wrong. run with the following command with the following command:
	tcsh 07_00_forLoop_checkProcessing_MNI.sh

08 run code 08_00_forLoop_checkProcessing_ORIG.sh this code runs creates some figure for the anatomy and EPI data for each participant (sub-01 to sub-26), calling the script '08_01_visualizeProcessed_ORIG.sh' for each participant. Results are stored in the folder 'derivatives/images_processing_afni_ORIG_denoised/'. These figures are from ORIG data, not warped into MNI. These figures are useful to quickly check whether some basic processing in terms of volreg and aligh went wrong. run with the following command with the following command:
	tcsh 08_00_forLoop_checkProcessing_ORIG.sh

09 get rois from atlases (Glasser and Princeton visual). Brings these atlases into the individual participant Freesurfer folder
	sh 09_00_getRoisInIndividualSpace.sh

10 get rois from atlases (Glasser and Princeton visual). Brings these atlases into the MNI space in volumetric coordinates (not surface based) MNI folder: /media/alessiofracasso/DATADRIVE1/Flanker/derivatives/surfaceAtlases/suma_MNI152_2009_princetonAtlas/
	sh 10_00_getRoisInMNISpace.sh

11 11_00_copyIndividualRois_ORIG_DENOISED.sh. Copies ROIs in ORIG space from Freesurfer individual folder, to individual processing folders: derivatives/processing_afni_denoised. ROI data in interpolated with nearest neighbour to the stat bucket in the processing folder 

12 12_00_copyIndividualRois_ORIG_NO_DENOISED.sh. Copies ROIs in ORIG space from Freesurfer individual folder, to individual processing folders: derivatives/processing_afni_no_denoised. ROI data in interpolated with nearest neighbour to the stat bucket in the processing folder 

13 13_00_copyIndividualRois_MNI_DENOISED.sh. Copies ROIs in MNI space from Freesurfer MNI folder, to individual processing folders: derivatives/processing_afni_MNI_denoised. ROI data in interpolated with nearest neighbour to the stat bucket in the processing folder

14 14_00_copyIndividualRois_MNI_no_DENOISED.sh. Copies ROIs in MNI space from Freesurfer MNI folder, to individual processing folders: derivatives/processing_afni_MNI_no_denoised. ROI data in interpolated with nearest neighbour to the stat bucket in the processing folder

15_volumetricStats_Glasser.Rextract volumetric stats from Glasser and save output

16_analysis_Glasser.R analysis Glasser

17_volumetricStats_Wang.R volumetric stats from Wang and save output

18_analysis_Wang.R Wang analysis

##############################


03 run code stat_3dTtest_command.sh with the following command:
	tcsh stat_3dTtest_command.sh
this code runs some stats on the incongruent - congruent comparison across participants using 3dtstat++ and 3dMEMA

04 runFreesurferSuma.R
	runs freesurfer and suma on the anatomy of each participant
	
05 run code forLoop_checkSegmentations.sh with the following command:
	tcsh forLoop_checkSegmentations.sh

this code runs creates some figure for the anatomy segmentation for each participant (sub-01 to sub-26), calling the script 'visualizeSegmentations.sh' for each participant.
results are stored in the folder 'derivatives/images_processing_images_processing_Freesurfer/'. These figures are useful to quickly check whether some basic segmentation gave problems.

