#!/bin/bash

3dDeconvolve -input pb04.$subj.r*.scale+orig.HEAD                        \
    -censor motion_${subj}_censor.1D                                     \
    -ortvec mot_demean.r01.1D mot_demean_r01                             \
    -ortvec mot_demean.r02.1D mot_demean_r02                             \
    -polort 2 -float                                                     \
    -basis_normall 1.0                                                   \
    -num_stimts 1                                                        \
    -stim_times_AM1 1 stimuli/RT_correct.1D 'dmUBLOCK'                   \
    -stim_label 1 RT_correct                                             \
    -jobs 4                                                              \
    -GOFORIT 10                                                          \
    -fout -tout -x1D X.xmat.1D -xjpeg X.jpg                              \
    -x1D_uncensored X.nocensor.xmat.1D                                   \
    -fitts fitts.$subj                                                   \
    -errts errts.${subj}                                                 \
    -bucket stats.$subj



# run the regression analysis
3dDeconvolve -input pb04.$subj.r*.scale+orig.HEAD                        \
    -censor motion_${subj}_censor.1D                                     \
    -ortvec mot_demean.r01.1D mot_demean_r01                             \
    -ortvec mot_demean.r02.1D mot_demean_r02                             \
    -polort 2 -float                                                     \
    -basis_normall 1.0                                                   \
    -num_stimts 1                                                        \
    -stim_times_IM 1 stimuli/RT_correctWrong_all_for_IM.1D 'dmUBLOCK'    \
    -stim_label 1 stim01                                                 \
    -jobs 4                                                              \
    -GOFORIT 10                                                          \
    -fout -tout -x1D X.xmat.1D -xjpeg X.jpg                              \
    -x1D_uncensored X.nocensor.xmat.1D                                   \
    -fitts fitts.$subj                                                   \
    -errts errts.${subj}                                                 \
    -bucket stats.$subj


# ------------------------------
# run the regression analysis
3dDeconvolve -input pb04.$subj.r*.scale+orig.HEAD                        \
    -censor motion_${subj}_censor.1D                                     \
    -ortvec mot_demean.r01.1D mot_demean_r01                             \
    -ortvec mot_demean.r02.1D mot_demean_r02                             \
    -polort 2 -float                                                     \
    -num_stimts 2                                                        \
    -stim_times 1 stimuli/congruent.1D 'BLOCK(1,1)'                      \
    -stim_label 1 congruent                                              \
    -stim_times 2 stimuli/incongruent.1D 'BLOCK(1,1)'                    \
    -stim_label 2 incongruent                                            \
    -jobs 8                                                              \
    -gltsym 'SYM: incongruent -congruent'                                \
    -glt_label 1 incongruent-congruent                                   \
    -gltsym 'SYM: congruent -incongruent'                                \
    -glt_label 2 congruent-incongruent                                   \
    -fout -tout -x1D X.xmat.1D -xjpeg X.jpg                              \
    -x1D_uncensored X.nocensor.xmat.1D                                   \
    -fitts fitts.$subj                                                   \
    -errts errts.${subj}                                                 \
    -bucket stats.$subj


