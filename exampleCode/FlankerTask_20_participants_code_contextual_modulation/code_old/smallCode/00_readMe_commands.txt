unzip the folders suma_MNI152_2009_princetonAtlas the derivatives/surfaceAtlas 
the folder suma_MNI152_2009_princetonAtlas contains the MNI surfaces and the atlas information, you can visualize the MNI volume and the MNI surfaces using the command 'afni -niml & suma -spec std.141.MNI152_2009_both.spec -sv MNI152_2009_SurfVol.nii' from a terminal in suma_MNI152_2009_princetonAtlas. Once the surfaces are open, you can right click on the surface and press 't' on the keyboard to link the surfaces with the afni volume. Then you can load the Wang 2015 atlas, follow the video I sent you.


03 run code 03_00_forLoop_MNI_no_denoised.sh with the command below. This code runs the main processing across all participants (sub-01 to sub-26). Results are stored in the folder 'derivatives/processing_afni_MNI_no_denoised/'. This code calls the code: 03_01_sub_xx_afni_proc_MNI_no_denoised.sh.
to run, from the main folder:
	# tcsh 03_00_forLoop_MNI_no_denoised.sh

10 get rois from atlases (Glasser and Princeton visual). Brings these atlases into the MNI space in volumetric coordinates (not surface based) MNI folder: Flanker/derivatives/surfaceAtlases/suma_MNI152_2009_princetonAtlas/
remember to adjust the folders set in the file according to your path before running it
	sh 10_00_getRoisInMNISpace.sh

14 14_00_copyIndividualRois_MNI_no_DENOISED.sh. Copies ROIs in MNI space from Freesurfer MNI folder, to individual processing folders: derivatives/processing_afni_MNI_no_denoised. ROI data in interpolated with nearest neighbour to the stat bucket in the processing folder
remember to adjust the folders set in the file according to your path before running it

17_volumetricStats_Wang.R volumetric stats from Wang and save output. Output is saved in the folder derivatives/resultsWang, which needs to be created manually. To run the code you can type 
Rscript 17_volumetricStats_Wang.R from the terminal open in your code folder. Remember to change the folder names according to your path in the file before running it

18_analysis_Wang.R Wang analysis

##############################

