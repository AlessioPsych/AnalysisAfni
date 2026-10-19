#!/usr/bin/env tcsh

# set subject and group identifiers
set subj  = $argv[1]
set main_folder = $argv[2]
set input_folder = $argv[3]
set gname = Flanker

echo "subj function loop, main folder: '$main_folder'"
echo "subj function loop, input folder: '$input_folder'"


# set data directories
set top_dir = $input_folder/{$subj}
set anat_dir  = $top_dir/anat
set epi_dir   = $top_dir/func
set stim_dir  = $top_dir/func

# run afni_proc.py to create a single subject processing script
#afni_proc.py -subj_id $subj                                      \
#        -script proc.$subj -scr_overwrite                        \
#        -blocks tshift 					\
#        -copy_anat $anat_dir/{$subj}_T1w.nii.gz                   \
#        -dsets                                                   \
#            $epi_dir/{$subj}_task-flanker_run-1_bold.nii.gz       \
#            $epi_dir/{$subj}_task-flanker_run-2_bold.nii.gz \
#        -no_epi_review -html_review_style none -execute
   


# run afni_proc.py to create a single subject processing script
afni_proc.py -subj_id $subj                                      \
        -script proc.$subj -scr_overwrite                        \
        -blocks tshift align volreg blur mask scale regress \
        #-blocks tshift volreg blur mask scale regress \
        -copy_anat $anat_dir/{$subj}_T1w.nii.gz                  \
        -dsets                                                   \
            $epi_dir/{$subj}_task-flanker_run-1_bold_denoised.nii.gz       \
            $epi_dir/{$subj}_task-flanker_run-2_bold_denoised.nii.gz       \
        -tcat_remove_first_trs 0                                 \
        -align_opts_aea -cost lpc -giant_move                    \
        -volreg_align_to MIN_OUTLIER                             \
        -volreg_align_e2a                                        \
        -blur_size 3                                             \
        -regress_stim_times                                      \
            $stim_dir/congruent.1D                               \
            $stim_dir/incongruent.1D                             \
            $stim_dir/RT_all.1D                                  \
        -regress_stim_labels                                     \
            congruent incongruent rt 				 \
        #-regress_stim_types times times AM2                      \
        #-regress_basis_multi 'BLOCK(1,1) BLOCK(1,1) dmBLOCK'     \
        -regress_stim_types times times AM1                      \
        -regress_basis_multi 'BLOCK(0.1,1) BLOCK(0.1,1) dmUBLOCK'     \
        -regress_censor_motion 0.3                               \
        -regress_motion_per_run                                  \
        -regress_opts_3dD                                        \
            -jobs 4 -GOFORIT 10                                  \
            -gltsym 'SYM: incongruent -congruent' -glt_label 1   \
        incongruent-congruent                                    \
	-gltsym 'SYM: congruent -incongruent' -glt_label 2       \
        congruent-incongruent                                    \
        -regress_reml_exec -regress_opts_reml -GOFORIT           \
        -regress_make_ideal_sum sum_ideal.1D                     \
        -regress_est_blur_epits                                  \
        -regress_est_blur_errts                                  \
        -regress_run_clustsim no -no_epi_review -html_review_style none -remove_preproc_files -execute


