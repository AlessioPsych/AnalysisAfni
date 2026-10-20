#!/bin/tcsh -xef

set main_folder = '/home/fracasso/Data/openNeuro/ds001751_FlankerContextualModulation'
set code_folder = $main_folder/code
set input_folder = $main_folder/derivatives/mrtrix3
set output_folder = $main_folder/derivatives/processing_afni_denoised_incongruent_congruent_tent
echo $main_folder
echo $code_folder
echo $input_folder
echo $output_folder

cd $input_folder

# creates output folder

#if ( -d $output_folder ) then
#    echo "Directory '$output_folder' exists. Removing and recreating..."
#    rm -rf $output_folder
#else
#    echo "Directory '$output_folder' does not exist. Creating..."
#endif

#mkdir $output_folder

# creates subjList.txt file
if (! -f subjList.txt) then
    ls | grep ^sub- > subjList.txt
endif

# runs code across participants
foreach i (`cat subjList.txt`)
#set i = 'sub-01'
	
	set subj = $i
	cd $input_folder
	cd {$i}/func
	pwd

	cd $output_folder
	cd {$i}.results
	pwd

	ls *_RT_correct_basis_1*
	rm *_RT_correct_basis_1*

	# -------------------------------------------
	# run the regression analysis, RT tent approach
	# -------------------------------------------

	3dDeconvolve -input pb04.$subj.r*.scale+orig.HEAD                        \
	    -censor motion_${subj}_censor.1D                                     \
	    -ortvec mot_demean.r01.1D mot_demean_r01                             \
    	-ortvec mot_demean.r02.1D mot_demean_r02                             \
    	-polort 9 -float                                                     \
    	-num_stimts 4                                                        \
    	-basis_normall 1.0                                                   \
    	-stim_times 1 stimuli/cue.1D 'BLOCK(1,1)'                            \
    	-stim_label 1 cue                                                    \
    	-stim_times 2 stimuli/blank_MC.1D 'TENT(0,15,12)'                       \
    	-stim_label 2 blank_MC                                                  \
    	-stim_times 3 stimuli/blank_MI.1D 'TENT(0,15,12)'                       \
    	-stim_label 3 blank_MI                                                  \
	-stim_times_AM1 4 $input_folder/{$i}/func/RT_correct.1D 'dmUBLOCK'   \
	-stim_label 4 RT_correct                                             \
    	-jobs 8                                                              \
    	-fout -tout -x1D X.xmat_RT_correct_basis_1.1D -xjpeg X_RT_correct_basis_1.jpg                              \
    	-x1D_uncensored X.nocensor.xmat_RT_correct_basis_1.1D                                   \
    	-fitts fitts_RT_correct_basis_1.$subj                                                   \
    	-errts errts_RT_correct_basis_1.${subj}                                                 \
    	-bucket stats_RT_correct_basis_1.$subj


	ls *_RT_correct_basis_no*
	rm *_RT_correct_basis_no*

	# -------------------------------------------
	# run the regression analysis, RT tent approach (no)
	# -------------------------------------------

	3dDeconvolve -input pb04.$subj.r*.scale+orig.HEAD                        \
	    -censor motion_${subj}_censor.1D                                     \
	    -ortvec mot_demean.r01.1D mot_demean_r01                             \
    	-ortvec mot_demean.r02.1D mot_demean_r02                             \
    	-polort 9 -float                                                     \
    	-num_stimts 4                                                        \
    	-stim_times 1 stimuli/cue.1D 'BLOCK(1,1)'                            \
    	-stim_label 1 cue                                                    \
    	-stim_times 2 stimuli/blank_MC.1D 'TENT(0,15,12)'                       \
    	-stim_label 2 blank_MC                                                  \
    	-stim_times 3 stimuli/blank_MI.1D 'TENT(0,15,12)'                       \
    	-stim_label 3 blank_MI                                                  \
	-stim_times_AM1 4 $input_folder/{$i}/func/RT_correct.1D 'dmUBLOCK'   \
	-stim_label 4 RT_correct                                             \
    	-jobs 8                                                              \
    	-fout -tout -x1D X.xmat_RT_correct_basis_no.1D -xjpeg X_RT_correct_basis_no.jpg                              \
    	-x1D_uncensored X.nocensor.xmat_RT_correct_basis_no.1D                                   \
    	-fitts fitts_RT_correct_basis_no.$subj                                                   \
    	-errts errts_RT_correct_basis_no.${subj}                                                 \
    	-bucket stats_RT_correct_basis_no.$subj








	
#	ls *_RT_correct_basis_1*
#	rm *_RT_correct_basis_1*
	# RT_correct -basis_normall 1.0
#	3dDeconvolve -input pb04.$subj.r*.scale+orig.HEAD                        \
#	    -censor motion_${subj}_censor.1D                                     \
#	    -ortvec mot_demean.r01.1D mot_demean_r01                             \
#	    -ortvec mot_demean.r02.1D mot_demean_r02                             \
#	    -polort 2 -float                                                     \
#	    -basis_normall 1.0                                                   \
#	    -num_stimts 1                                                        \
#	    -stim_times_AM1 1 $input_folder/{$i}/func/RT_correct.1D 'dmUBLOCK'   \
#	    -stim_label 1 RT_correct                                             \
#	    -jobs 8                                                              \
#	    -GOFORIT 10                                                          \
#	    -fout -tout -x1D X_RT_correct_basis_1.xmat.1D -xjpeg X_RT_correct_basis_1.jpg          \
#	    -x1D_uncensored X_RT_correct_basis_1.nocensor.xmat.1D                                  \
#	    -fitts fitts_RT_correct_basis_1.$subj                                                  \
#	    -errts errts_RT_correct_basis_1.${subj}                                                \
#	    -bucket stats_RT_correct_basis_1.$subj

#	ls *_RT_correct_basis_no*
#	rm *_RT_correct_basis_no*
	# RT_correct NO -basis_normall 1.0
#	3dDeconvolve -input pb04.$subj.r*.scale+orig.HEAD                        \
#	    -censor motion_${subj}_censor.1D                                     \
#	    -ortvec mot_demean.r01.1D mot_demean_r01                             \
#	    -ortvec mot_demean.r02.1D mot_demean_r02                             \
#	    -polort 2 -float                                                     \
#	    -num_stimts 1                                                        \
#	    -stim_times_AM1 1 $input_folder/{$i}/func/RT_correct.1D 'dmUBLOCK'   \
#	    -stim_label 1 RT_correct                                             \
#	    -jobs 8                                                              \
#	    -GOFORIT 10                                                          \
#	    -fout -tout -x1D X_RT_correct_basis_no.xmat.1D -xjpeg X_RT_correct_basis_no.jpg         \
#	    -x1D_uncensored X_RT_correct_basis_no.nocensor.xmat.1D                                  \
#	    -fitts fitts_RT_correct_basis_no.$subj                                                  \
#	    -errts errts_RT_correct_basis_no.${subj}                                                \
#	    -bucket stats_RT_correct_basis_no.$subj


    echo ................
    echo ................
    echo ................
    
end




