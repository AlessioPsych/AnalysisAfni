#!/usr/bin/env tcsh

# set subject and group identifiers
set subj  = $argv[1]
set main_folder = $argv[2]
set input_folder = $argv[3]
set anat_folder = $argv[4]
set gname = Flanker

echo "subj function loop, main folder: '$main_folder'"
echo "subj function loop, input folder: '$input_folder'"


# set data directories
set top_dir = $input_folder/{$subj}
set anat_dir  = $anat_folder/{$subj}/Freesurfer_result/SUMA
set epi_dir   = $top_dir/func
set stim_dir  = $top_dir/func

#if ( "$cond1" == "something" || "$cond1" == "something_else" ) then
    # Do something
#else
    # Do something else
#endif

# run afni_proc.py to create a single subject processing script

#if ( $subj != "subj-03" && $subj != "subj-16" ) then

#afni_proc.py -subj_id $subj                                      \
#        -script proc.$subj -scr_overwrite                        \
#        -blocks align volreg \
#        -copy_anat $anat_dir/brainAFNI+orig                  \
#        -dsets                                                   \
#            $epi_dir/{$subj}_task-flanker_run-1_bold_denoised_tSliceCorrect.nii.gz       \
#            $epi_dir/{$subj}_task-flanker_run-2_bold_denoised_tSliceCorrect.nii.gz       \
#        -tcat_remove_first_trs 0                                 \
#        -align_opts_aea -cost lpc -ginormous_move -epi2anat -deoblique off -anat_has_skull no \
#        -volreg_align_to MIN_OUTLIER                             \
#        -volreg_align_e2a                                        \
#        -no_epi_review -html_review_style none -execute
#else
#       #-Allineate_opts -twopass -weight_frac 1.0 -maxrot 8 -maxshf 15 -warp aff -source_automask+4
 
#afni_proc.py -subj_id $subj                                      \
#        -script proc.$subj -scr_overwrite                        \
#        -blocks align volreg \
#        -copy_anat $anat_dir/brainAFNI+orig                  \
#        -dsets                                                   \
#            $epi_dir/{$subj}_task-flanker_run-1_bold_denoised_tSliceCorrect.nii.gz       \
#            $epi_dir/{$subj}_task-flanker_run-2_bold_denoised_tSliceCorrect.nii.gz       \
#        -tcat_remove_first_trs 0                                 \
#        -align_opts_aea -cost lpc -giant_move -deoblique on \
#        -volreg_align_to MIN_OUTLIER                             \
#        -volreg_align_e2a                                        \
#        -no_epi_review -html_review_style none -execute
#endif

   


# run afni_proc.py to create a single subject processing script
afni_proc.py -subj_id $subj                                      \
        -script proc.$subj -scr_overwrite                        \
        -blocks tshift align volreg blur mask scale regress \
        -copy_anat $anat_dir/brainAFNI+orig                  \
        -dsets                                                   \
            $epi_dir/{$subj}_task-flanker_run-1_bold_denoised_tSliceCorrect_deob.nii.gz       \
            $epi_dir/{$subj}_task-flanker_run-2_bold_denoised_tSliceCorrect_deob.nii.gz       \
        -tcat_remove_first_trs 0                                 \
        -align_opts_aea -cost lpc -ginormous_move \
        -volreg_align_to MIN_OUTLIER                             \
        -volreg_align_e2a                                        \
        -blur_size 2                                             \
        -regress_stim_times                                      \
            $stim_dir/cue.1D 					 \
            $stim_dir/blank_MC.1D 					 \
            $stim_dir/blank_MI.1D 					 \
            $stim_dir/congruent_MC.1D                               \
            $stim_dir/incongruent_MC.1D                             \
            $stim_dir/congruent_MI.1D                               \
            $stim_dir/incongruent_MI.1D                             \
        -regress_stim_labels                                     \
            cue blank_MC blank_MI congruent_MC incongruent_MC congruent_MI incongruent_MI                      \
        -regress_basis_multi 'BLOCK(1,1)' 'TENT(0,15,12)' 'TENT(0,15,12)' 'BLOCK(1,1)' 'BLOCK(1,1)' 'BLOCK(1,1)' 'BLOCK(1,1)'    \
        -regress_censor_motion 0.3                               \
        -regress_motion_per_run                                  \
        -regress_opts_3dD                                        \
            -jobs 10                                              \
            -gltsym 'SYM: +incongruent_MC +incongruent_MI -congruent_MC -congruent_MI' -glt_label 1   \
        incongruent-congruent                                    \
	-gltsym 'SYM: +congruent_MC +congruent_MI -incongruent_MC -incongruent_MI' -glt_label 2       \
        congruent-incongruent                                    \
        -regress_reml_exec                                       \
        -regress_make_ideal_sum sum_ideal.1D                     \
        -regress_est_blur_epits                                  \
        -regress_est_blur_errts                                  \
        -regress_run_clustsim no -no_epi_review -html_review_style none -execute


