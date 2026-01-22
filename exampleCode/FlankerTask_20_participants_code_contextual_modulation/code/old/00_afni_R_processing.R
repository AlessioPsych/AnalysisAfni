subjId <- 'sub-08'
mainDir <- sprintf( '/media/alessiofracasso/DATADRIVE1/Flanker' )
anatDir <- sprintf( '%s/anat', subjId )
epiDir <-  sprintf( '%s/func', subjId )
stimDir <-  sprintf( '%s/func', subjId )

# run afni_proc.py to create a single subject processing script

instr <- paste( sprintf('afni_proc.py -subj_id %s', subjId),
                sprintf('-script proc.%s -scr_overwrite',subjId),
                sprintf('-copy_anat %s/%s_T1w.nii.gz', anatDir, subjId ),
                '-blocks tshift align tlrc volreg mask blur scale',
                '-dsets',
                sprintf('%s/%s_task-flanker_run-1_bold.nii.gz', epiDir, subjId),
                sprintf('%s/%s_task-flanker_run-2_bold.nii.gz', epiDir, subjId),
                '-tcat_remove_first_trs 0',
                '-align_unifize_epi local',                          	
                '-align_opts_aea -cost lpc -giant_move',
                '-tlrc_base MNI152_2009_template.nii.gz',
                '-tlrc_opts_at -init_xform AUTO_CENTER',
                '-volreg_align_e2a',
                '-volreg_tlrc_warp',     
                '-volreg_compute_tsnr yes',
                '-mask_epi_anat yes',
                '-blur_size 4.0', 
                '-no_epi_review', 
                '-html_review_style none',
                '-execute' )
setwd( mainDir  )
print( getwd() )
print( instr )
system( instr )

