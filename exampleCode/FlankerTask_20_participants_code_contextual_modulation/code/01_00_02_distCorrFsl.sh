#!/bin/bash

# Instructions to run:
# sudo docker run -v /mnt/disk01/ds001751_FlankerTask_context:/mrtrixDataFolder --rm -it alessiodock/fsl_60710:01

# sudo docker run -it fslimage:Dockerfile
# sudo docker run -e DISPLAY=$DISPLAY  -v /tmp/.X11-unix:/tmp/.X11-unix -v /home/fracasso/Data/openNeuro/ds001751:/mrtrixDataFolder --rm --user $(id -u):$(id -g) -it fslimage:Dockerfile

# cd /mrtrixDataFolder
# sh 01_00_02_distCorrFsl.sh

basedir="/mrtrixDataFolder"
maindir="$basedir/derivatives/mrtrix3"
#targetdir="derivatives/mrtrix3"

#apt-get install bc

echo "$maindir"

cd "$maindir"
echo "Current folder: $PWD"

# print list of participants

# to run properly on all participants
#\ls -d sub-* > subdirs.txt

# to test only:
\ls -d sub-01 sub-02 > subdirs.txt
chmod 777 subdirs.txt
#dir="sub-01"

#if [ -d $maindir/$targetdir ]; then
#	echo "removing folder $maindir/$targetdir" 
#	rm -R $maindir/$targetdir
#fi
#mkdir $maindir/$targetdir

while IFS= read -r dir; do

	cd "$maindir/$dir/fmap"
	echo "Processing directory: $dir"
	echo
	echo "folder: $PWD"
	echo
	
	
	raw_dir=$maindir/$dir/fmap
	output_dir=$maindir/$dir/fmap
	anat_dir=$maindir/$dir/anat
	func_dir=$maindir/$dir/func

	# Reorient fmap magnitude and phase image
	
	echo "if [ -f fmap_mag.nii.gz ]; then"
	echo "rm fmap_mag.nii.gz"
  	echo "fslreorient2std ${dir}_magnitude2.nii.gz fmap_mag.nii.gz"
	if [ -f fmap_mag.nii.gz ]; then
		rm fmap_mag.nii.gz
	fi
	fslreorient2std ${dir}_magnitude2.nii.gz fmap_mag.nii.gz
	chmod 777 fmap_mag.nii.gz

	echo "if [ -f fmap_phase.nii.gz ]; then"
	echo "rm fmap_phase.nii.gz"
  	echo "fslreorient2std ${dir}_phasediff.nii.gz fmap_phase.nii.gz"
	if [ -f fmap_phase.nii.gz ]; then
		rm fmap_phase.nii.gz
	fi
	fslreorient2std ${dir}_phasediff.nii.gz fmap_phase.nii.gz
	chmod 777 fmap_phase.nii.gz

  	# Run bias correction on fmap magnitude image
	if [ -d fmap.anat ]; then
		rm -R fmap.anat/
	fi
	echo "Running fsl_anat on fmap"
  	fsl_anat --strongbias  --nocrop  --noreg  --nosubcortseg  --noseg -i fmap_mag.nii.gz -o fmap
	chmod 777 -R fmap.anat/
	
	# Bet fmap magnitude image
	if [ -f fmap.anat/fmap_mag_brain_bet.nii.gz ]; then
		rm fmap.anat/fmap_mag_brain_bet.nii.gz
	fi
	if [ -f fmap.anat/fmap_mag_brain_bet.nii.gz ]; then
		rm fmap.anat/fmap_mag_brain_bet.nii
	fi
	bet fmap.anat/T1_biascorr.nii.gz fmap.anat/fmap_mag_brain_bet.nii.gz -R
	chmod 777 fmap.anat/fmap_mag_brain_bet.nii.gz

	if [ -f fmap.anat/map_mag_brain_ero1.nii.gz ]; then
		rm fmap.anat/map_mag_brain_ero1.nii.gz
	fi
	if [ -f fmap.anat/map_mag_brain_ero2.nii.gz ]; then
		rm fmap.anat/map_mag_brain_ero2.nii.gz
	fi	
	fslmaths fmap.anat/fmap_mag_brain_bet.nii.gz -ero fmap.anat/fmap_mag_brain_ero1.nii.gz
	chmod 777 fmap.anat/fmap_mag_brain_ero1.nii.gz
	fslmaths fmap.anat/fmap_mag_brain_ero1.nii.gz -ero fmap.anat/fmap_mag_brain_ero2.nii.gz
	chmod 777 fmap.anat/fmap_mag_brain_ero2.nii.gz

	chosen_brain_extraction=fmap.anat/fmap_mag_brain_ero1.nii.gz
	
	# Copy the skull-stripped and non-stripped magnitude image and rename them
	if [ -f fmap_mag_brain.nii.gz ]; then
		rm fmap_mag_brain.nii.gz
	fi	
	if [ -f fmap_mag.nii.gz ]; then
		rm fmap_mag.nii.gz
	fi	
	fslmaths ${chosen_brain_extraction} ${output_dir}/fmap_mag_brain.nii.gz
	fslmaths fmap.anat/T1_biascorr ${output_dir}/fmap_mag.nii.gz
	chmod 777 fmap_mag_brain.nii.gz
	chmod 777 fmap_mag.nii.gz
	
	# rescale phase data
	if [ -f fmap_phase_scaled.nii.gz ]; then
		rm fmap_phase_scaled.nii.gz
	fi	
	infile=fmap_phase.nii.gz
	outfile=fmap_phase_scaled.nii.gz

	# get min and max
	read a b <<< $(fslstats "$infile" -R)

	# compute scale factor 2π/(b-a)
	scale=$(awk -v a="$a" -v b="$b" 'BEGIN {print 6.283185307179586/(b-a)}')

	# rescale: (x - a) * scale
	fslmaths "$infile" -sub "$a" -mul "$scale" "$outfile"
	chmod 777 fmap_phase_scaled.nii.gz

	echo "Input range: $a  $b"
	echo "Scale factor: $scale"
	echo "Output written to $outfile"
	echo "Output file range:"
	fslstats "$outfile" -R
	
	# Generate fieldmap, I have bc (basic computations) installed, so it will scale the phase to the correct range (0 - 6.28)
	if [ -f fmap_rads.nii.gz ]; then
		rm fmap_rads.nii.gz
	fi	
	fsl_prepare_fieldmap SIEMENS ${output_dir}/fmap_phase.nii.gz ${output_dir}/fmap_mag_brain.nii.gz ${output_dir}/fmap_rads.nii.gz 2.46
	chmod 777 fmap_rads.nii.gz
	
	## clean up
	#rm fmap_*
	#rm -R fmap.anat/

	# prepare T1 image
	if [ -f ${dir}_T1w.nii.gz ]; then
		rm ${dir}_T1w.nii.gz
	fi	
	cp $anat_dir/${dir}_T1w.nii.gz ${output_dir}/${dir}_T1w.nii.gz
	chmod 777 ${dir}_T1w.nii.gz

	if [ -f fmap_T1w_brain.nii.gz ]; then
		rm fmap_T1w_brain.nii.gz
	fi		
	bet ${dir}_T1w.nii.gz fmap_T1w_brain.nii.gz
	chmod 777 fmap_T1w_brain.nii.gz
	
	rm fmap_T1w_brain_*
	fast -B -I 10 -l 10 fmap_T1w_brain.nii.gz
	chmod 777 fmap_T1w_brain*
	
	rm fmap_T1w_wmseg.nii.gz
	fslmaths fmap_T1w_brain_pve_2.nii.gz -thr 0.5 -bin fmap_T1_wmseg.nii.gz
	chmod 777 fmap_T1_wmseg.nii.gz

	cp $func_dir/${dir}_task-flanker_run-1_bold_denoised_tSliceCorrect.nii.gz ${output_dir}/run-1_bold.nii.gz
	chmod 777 run-1_bold.nii.gz
	cp $func_dir/${dir}_task-flanker_run-2_bold_denoised_tSliceCorrect.nii.gz ${output_dir}/run-2_bold.nii.gz
	chmod 777 run-2_bold.nii.gz
	
	# provide reasonable coregistration
	epidata=run-1_bold.nii.gz
	flirt -in $epidata -ref fmap_T1w_brain.nii.gz -out func_in_struct_start.nii.gz -omat func_to_struct.mat -dof 6 -cost normmi
	chmod 777 func_in_struct_start*
	chmod 777 func_to_struct*
	
	# apply reasonable start to rads, mag, and stripped mag:
	flirt -in fmap_rads.nii.gz -ref fmap_T1w_brain.nii.gz -out fmap_rads_coreg.nii.gz -init func_to_struct.mat -applyxfm
	chmod 777 fmap_rads_coreg.nii.gz

	flirt -in fmap_mag.nii.gz -ref fmap_T1w_brain.nii.gz -out fmap_mag_coreg.nii.gz -init func_to_struct.mat -applyxfm
	chmod 777 fmap_mag_coreg.nii.gz

	flirt -in fmap_mag_brain.nii.gz -ref fmap_T1w_brain.nii.gz -out fmap_mag_brain_coreg.nii.gz -init func_to_struct.mat -applyxfm
	chmod 777 fmap_mag_brain_coreg.nii.gz
	
	# BOLD run 01
	epi_reg --echospacing=0.00059 --wmseg=fmap_T1_wmseg.nii.gz --fmap=fmap_rads_coreg.nii.gz --fmapmag=fmap_mag_coreg.nii.gz --fmapmagbrain=fmap_mag_brain_coreg.nii.gz --pedir=-y --epi=func_in_struct_start.nii.gz --t1=${dir}_T1w.nii.gz --t1brain=fmap_T1w_brain.nii.gz --out=func2struct_negy
	chmod 777 func2struct_negy*

	epi_reg --echospacing=0.00059 --wmseg=fmap_T1_wmseg.nii.gz --fmap=fmap_rads_coreg.nii.gz --fmapmag=fmap_mag_coreg.nii.gz --fmapmagbrain=fmap_mag_brain_coreg.nii.gz --pedir=y --epi=func_in_struct_start.nii.gz --t1=${dir}_T1w.nii.gz --t1brain=fmap_T1w_brain.nii.gz --out=func2struct_posy
	chmod 777 func2struct_posy*

	# apply basic coreg
	flirt -in run-1_bold.nii.gz -ref fmap_T1w_brain.nii.gz -out run-1_bold_coreg.nii.gz -init func_to_struct.mat -applyisoxfm 3
	chmod 777 run-1_bold_coreg.nii.gz

	# apply warp
	applywarp -i run-1_bold_coreg.nii.gz -r run-1_bold_coreg.nii.gz -o run-1_bold_coreg_warp.nii.gz -w func2struct_negy_warp.nii.gz 
	#--postmat=warpexample_func_bold01.mat
	chmod 777 run-1_bold_coreg_warp.nii.gz

	rm func_to_struct* 
	rm func2struct_posy*	
	rm func2struct_negy* 
	rm func_in_struct_start* 
	rm fmap_rads_coreg.nii.gz
	rm fmap_mag_brain_coreg.nii.gz
	rm fmap_mag_coreg.nii.gz


	epidata=run-2_bold.nii.gz
	flirt -in $epidata -ref fmap_T1w_brain.nii.gz -out func_in_struct_start.nii.gz -omat func_to_struct.mat -dof 6 -cost normmi
	chmod 777 func_in_struct_start*
	chmod 777 func_to_struct*

	flirt -in fmap_rads.nii.gz -ref fmap_T1w_brain.nii.gz -out fmap_rads_coreg.nii.gz -init func_to_struct.mat -applyxfm
	chmod 777 fmap_rads_coreg.nii.gz

	flirt -in fmap_mag.nii.gz -ref fmap_T1w_brain.nii.gz -out fmap_mag_coreg.nii.gz -init func_to_struct.mat -applyxfm
	chmod 777 fmap_mag_coreg.nii.gz

	flirt -in fmap_mag_brain.nii.gz -ref fmap_T1w_brain.nii.gz -out fmap_mag_brain_coreg.nii.gz -init func_to_struct.mat -applyxfm
	chmod 777 fmap_mag_brain_coreg.nii.gz
	
	# BOLD run 02
	epi_reg --echospacing=0.00059 --wmseg=fmap_T1_wmseg.nii.gz --fmap=fmap_rads_coreg.nii.gz --fmapmag=fmap_mag_coreg.nii.gz --fmapmagbrain=fmap_mag_brain_coreg.nii.gz --pedir=-y --epi=func_in_struct_start.nii.gz --t1=${dir}_T1w.nii.gz --t1brain=fmap_T1w_brain.nii.gz --out=func2struct_negy
	chmod 777 func2struct_negy*

	epi_reg --echospacing=0.00059 --wmseg=fmap_T1_wmseg.nii.gz --fmap=fmap_rads_coreg.nii.gz --fmapmag=fmap_mag_coreg.nii.gz --fmapmagbrain=fmap_mag_brain_coreg.nii.gz --pedir=y --epi=func_in_struct_start.nii.gz --t1=${dir}_T1w.nii.gz --t1brain=fmap_T1w_brain.nii.gz --out=func2struct_posy
	chmod 777 func2struct_posy*

	# apply basic coreg
	flirt -in run-2_bold.nii.gz -ref fmap_T1w_brain.nii.gz -out run-2_bold_coreg.nii.gz -init func_to_struct.mat -applyisoxfm 3
	chmod 777 run-2_bold_coreg.nii.gz

	# apply warp
	applywarp -i run-2_bold_coreg.nii.gz -r run-2_bold_coreg.nii.gz -o run-2_bold_coreg_warp.nii.gz -w func2struct_negy_warp.nii.gz 
	#--postmat=warpexample_func_bold01.mat
	chmod 777 run-2_bold_coreg_warp.nii.gz

	rm func_to_struct* 
	rm func2struct_posy*	
	rm func2struct_negy* 
	rm func_in_struct_start* 
	rm fmap_rads_coreg.nii.gz
	rm fmap_mag_brain_coreg.nii.gz
	rm fmap_mag_coreg.nii.gz
	
	#epi_reg --epi=run-1_bold.nii.gz --wmseg=fmap_T1_wmseg.nii.gz --t1=${dir}_T1w.nii.gz --t1brain=fmap_T1w_brain.nii.gz --out=func2struct_start
	#chmod 777 func2struct_start*
	
#	epi_reg --echospacing=0.00059 --wmseg=fmap_T1_wmseg.nii.gz --fmap=fmap_rads.nii.gz --fmapmag=fmap_mag.nii.gz --fmapmagbrain=fmap_mag_brain.nii.gz --pedir=-y --epi=run-1_bold.nii.gz --t1=${dir}_T1w.nii.gz --t1brain=fmap_T1w_brain.nii.gz --out=func2struct_negy
	
	echo "........."
	echo "........."
	echo "........."
          
done < subdirs.txt


