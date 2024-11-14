rm(list=ls()); gc()
for (nSubj in 1) {
  
  if ( nSubj==1 ) {
    mainDir <- '/analyse/Project0226/tests/nilsTest/data/20230622_SHI27'
    anatomyDir <- '/analyse/Project0226/tests/nilsTest/data/20230622_SHI27/ANATOMY/anat_test/SUMA'
  }
  
  setwd( mainDir )
  getwd()

  # model block design  
  # bring TS in anatomy  
  filenameIn <- 'ppp_modelOutput_blockDesign_DetrendedTs.nii.gz'
  filenameOut <- 'ppp_modelOutput_blockDesign_DetrendedTs_anatomy.nii.gz'
  if ( file.exists( sprintf('%s/%s', anatomyDir, filenameOut ) ) ) { system( sprintf('rm %s/%s', anatomyDir, filenameOut ) ) }
  instr <- paste( '3dAllineate',
                  sprintf('-prefix %s/%s', anatomyDir, filenameOut),
                  '-1Dmatrix_apply coregistration/coregMat.1D',
                  '-final lin',
                  sprintf('-input %s', filenameIn),
                  '-newgrid 2 2 2',
                  sprintf( '-master %s/anatCopy.nii.gz', anatomyDir ) )
  system( instr )
  
  # bring predictions in anatomy  
  filenameIn <- 'ppp_modelOutput_blockDesign_PredixtedTs.nii.gz'
  filenameOut <- 'ppp_modelOutput_blockDesign_PredixtedTs_anatomy.nii.gz'
  if ( file.exists( sprintf('%s/%s', anatomyDir, filenameOut ) ) ) { system( sprintf('rm %s/%s', anatomyDir, filenameOut ) ) }
  instr <- paste( '3dAllineate',
                  sprintf('-prefix %s/%s', anatomyDir, filenameOut),
                  '-1Dmatrix_apply coregistration/coregMat.1D',
                  '-final lin',
                  sprintf('-input %s', filenameIn),
                  '-newgrid 2 2 2',
                  sprintf( '-master %s/anatCopy.nii.gz', anatomyDir ) )
  system( instr )
  
  # bring params in anatomy  
  filenameIn <- 'ppp_modelOutput_blockDesign_params.nii.gz'
  filenameOut <- 'ppp_modelOutput_blockDesign_params_anatomy.nii.gz'
  if ( file.exists( sprintf('%s/%s', anatomyDir, filenameOut ) ) ) { system( sprintf('rm %s/%s', anatomyDir, filenameOut ) ) }
  instr <- paste( '3dAllineate',
                  sprintf('-prefix %s/%s', anatomyDir, filenameOut),
                  '-1Dmatrix_apply coregistration/coregMat.1D',
                  '-final lin',
                  sprintf('-input %s', filenameIn),
                  '-newgrid 2 2 2',
                  sprintf( '-master %s/anatCopy.nii.gz', anatomyDir ) )
  system( instr )
  
  # model block design and amplitude scaling
  # bring TS in anatomy  
  filenameIn <- 'ppp_modelOutput_ampScaling_DetrendedTs.nii.gz'
  filenameOut <- 'ppp_modelOutput_ampScaling_DetrendedTs_anatomy.nii.gz'
  if ( file.exists( sprintf('%s/%s', anatomyDir, filenameOut ) ) ) { system( sprintf('rm %s/%s', anatomyDir, filenameOut ) ) }
  instr <- paste( '3dAllineate',
                  sprintf('-prefix %s/%s', anatomyDir, filenameOut),
                  '-1Dmatrix_apply coregistration/coregMat.1D',
                  '-final lin',
                  sprintf('-input %s', filenameIn),
                  '-newgrid 2 2 2',
                  sprintf( '-master %s/anatCopy.nii.gz', anatomyDir ) )
  system( instr )
  
  # bring predictions in anatomy  
  filenameIn <- 'ppp_modelOutput_ampScaling_PredixtedTs.nii.gz'
  filenameOut <- 'ppp_modelOutput_ampScaling_PredixtedTs_anatomy.nii.gz'
  if ( file.exists( sprintf('%s/%s', anatomyDir, filenameOut ) ) ) { system( sprintf('rm %s/%s', anatomyDir, filenameOut ) ) }
  instr <- paste( '3dAllineate',
                  sprintf('-prefix %s/%s', anatomyDir, filenameOut),
                  '-1Dmatrix_apply coregistration/coregMat.1D',
                  '-final lin',
                  sprintf('-input %s', filenameIn),
                  '-newgrid 2 2 2',
                  sprintf( '-master %s/anatCopy.nii.gz', anatomyDir ) )
  system( instr )
  
  # bring params in anatomy  
  filenameIn <- 'ppp_modelOutput_ampScaling_params.nii.gz'
  filenameOut <- 'ppp_modelOutput_ampScaling_params_anatomy.nii.gz'
  if ( file.exists( sprintf('%s/%s', anatomyDir, filenameOut ) ) ) { system( sprintf('rm %s/%s', anatomyDir, filenameOut ) ) }
  instr <- paste( '3dAllineate',
                  sprintf('-prefix %s/%s', anatomyDir, filenameOut),
                  '-1Dmatrix_apply coregistration/coregMat.1D',
                  '-final lin',
                  sprintf('-input %s', filenameIn),
                  '-newgrid 2 2 2',
                  sprintf( '-master %s/anatCopy.nii.gz', anatomyDir ) )
  system( instr )
  
  setwd(anatomyDir)
  # get parameters volume into surface, simple block design, right hemisphere:
  filenameIn <- 'ppp_modelOutput_blockDesign_params_anatomy.nii.gz[5]'
  filenameOut <- 'ppp_modelOutput_blockDesign_params_anatomy.rh.niml.dset'
  if ( file.exists( sprintf('%s/%s', anatomyDir, filenameOut ) ) ) { system( sprintf('rm %s/%s', anatomyDir, filenameOut ) ) }
  instr <- paste( '3dVol2Surf',                 
                  '-spec std.141.anat_test_rh.spec',
                  '-surf_A std.141.rh.white.gii',
                  '-surf_B std.141.rh.pial.gii',
                  sprintf('-sv %s', filenameIn ),
                  sprintf('-grid_parent %s', filenameIn ),
                  '-map_func midpoint',
                  '-f_p1_fr 0',
                  '-f_pn_fr 0',
                  sprintf( '-out_niml %s', filenameOut ) )
  system( instr )
  
  # get parameters volume into surface, simple block design, left hemisphere:
  filenameIn <- 'ppp_modelOutput_blockDesign_params_anatomy.nii.gz[5]'
  filenameOut <- 'ppp_modelOutput_blockDesign_params_anatomy.lh.niml.dset'
  if ( file.exists( sprintf('%s/%s', anatomyDir, filenameOut ) ) ) { system( sprintf('rm %s/%s', anatomyDir, filenameOut ) ) }
  instr <- paste( '3dVol2Surf',                 
                  '-spec std.141.anat_test_lh.spec',
                  '-surf_A std.141.lh.white.gii',
                  '-surf_B std.141.lh.pial.gii',
                  sprintf('-sv %s', filenameIn ),
                  sprintf('-grid_parent %s', filenameIn ),
                  '-map_func midpoint',
                  '-f_p1_fr 0',
                  '-f_pn_fr 0',
                  sprintf( '-out_niml %s', filenameOut ) )
  system( instr )
  
  
}



# if ( file.exists( sprintf('%s/stats_concatenated_ts_masked_blurred.nii', anatomyDir ) ) ) { system( sprintf('rm %s/stats_concatenated_ts_masked_blurred.nii', anatomyDir ) ) }
# if ( file.exists( sprintf('%s/stat_mask_touch_min_rest_ts_concatenated.nii', anatomyDir ) ) ) { system( sprintf('rm %s/stat_mask_touch_min_rest_ts_concatenated.nii', anatomyDir ) ) }
# if ( file.exists( sprintf('%s/stat_active_min_rest_and_passive_min_rest_ts_concatenated.nii', anatomyDir ) ) ) { system( sprintf('rm %s/stat_active_min_rest_and_passive_min_rest_ts_concatenated.nii', anatomyDir ) ) }
# if ( file.exists( sprintf('%s/stat_mask_touch_min_rest_REML_ts_concatenated.nii', anatomyDir ) ) ) { system( sprintf('rm %s/stat_mask_touch_min_rest_REML_ts_concatenated.nii', anatomyDir ) ) }
# if ( file.exists( sprintf('%s/stats_concatenated_REML_ts_masked_blurred.nii', anatomyDir ) ) ) { system( sprintf('rm %s/stats_concatenated_REML_ts_masked_blurred.nii', anatomyDir ) ) }
# if ( file.exists( sprintf('%s/stat_active_min_rest_and_passive_min_rest_REML_ts_concatenated.nii', anatomyDir ) ) ) { system( sprintf('rm %s/stat_active_min_rest_and_passive_min_rest_REML_ts_concatenated.nii', anatomyDir ) ) }
# if ( file.exists( sprintf('%s/stat_active_plus_passive_min_touch_plus_instructions_ts_concatenated.nii', anatomyDir ) ) ) { system( sprintf('rm %s/stat_active_plus_passive_min_touch_plus_instructions_ts_concatenated.nii', anatomyDir ) ) }
# if ( file.exists( sprintf('%s/stat_active_plus_passive_min_touch_plus_instructions_REML_ts_concatenated.nii', anatomyDir ) ) ) { system( sprintf('rm %s/stat_active_plus_passive_min_touch_plus_instructions_REML_ts_concatenated.nii', anatomyDir ) ) }
# if ( file.exists( sprintf('%s/stat_active_min_passive_ts_concatenated.nii', anatomyDir ) ) ) { system( sprintf('rm %s/stat_active_min_passive_ts_concatenated.nii', anatomyDir ) ) }
# if ( file.exists( sprintf('%s/stat_active_min_passive_REML_ts_concatenated.nii', anatomyDir ) ) ) { system( sprintf('rm %s/stat_active_min_passive_REML_ts_concatenated.nii', anatomyDir ) ) }
# 
# 
# 
# 
# 
# # mask the stats on the anatomy only
# 
# 
# 
# # blur masked stats
# instr <- paste('3dmerge', 
#                '-1blur_fwhm 2',
#                '-doall',
#                sprintf('-prefix %s/stats_concatenated_ts_masked_blurred.nii', anatomyDir),
#                sprintf('%s/stats_concatenated_ts_masked.nii', anatomyDir) )
# system( instr )
# 
# # blur masked stats, REML
# instr <- paste('3dmerge', 
#                '-1blur_fwhm 2',
#                '-doall',
#                sprintf('-prefix %s/stats_concatenated_REML_ts_masked_blurred.nii', anatomyDir),
#                sprintf('%s/stats_concatenated_REML_ts_masked.nii', anatomyDir) )
# system( instr )
# 
# # conjunction analysis from concatenated ts #
# instr <- paste('3dcalc', # mask keeping only those voxels NOT responding to touch-rest,nor to instructions-rest, within the anatomy
#                sprintf('-a %s/anatCopy_resampled_to_stats.nii.gz', anatomyDir),
#                sprintf('-b %s/stats_concatenated_ts_masked_blurred.nii[23]', anatomyDir),
#                sprintf('-c %s/stats_concatenated_ts_masked_blurred.nii[26]', anatomyDir),
#                sprintf('-expr \u0027step(a)*not( step(abs(b)-4) )*not( step(abs(c)-4) )\u0027'),
#                sprintf('-prefix %s/stat_mask_touch_min_rest_ts_concatenated.nii', anatomyDir),
#                '-datum FLOAT' )
# system( instr )
# 
# # conjunction analysis from concatenated ts # REML
# instr <- paste('3dcalc', # mask keeping only those voxels NOT responding to touch-rest,nor to instructions-rest, within the anatomy
#                sprintf('-a %s/anatCopy_resampled_to_stats.nii.gz', anatomyDir),
#                sprintf('-b %s/stats_concatenated_REML_ts_masked_blurred.nii[23]', anatomyDir), # 'touch' tmap
#                sprintf('-c %s/stats_concatenated_REML_ts_masked_blurred.nii[26]', anatomyDir), # 'instructions' tmap
#                sprintf('-expr \u0027step(a)*not( step(abs(b)-4) )*not( step(abs(c)-4) )\u0027'),
#                sprintf('-prefix %s/stat_mask_touch_min_rest_REML_ts_concatenated.nii', anatomyDir),
#                '-datum FLOAT' )
# system( instr )
# 
# instr <- paste('3dcalc',
#                sprintf('-a %s/stats_concatenated_ts_masked_blurred.nii[17]', anatomyDir), # 'active' tmap
#                sprintf('-b %s/stats_concatenated_ts_masked_blurred.nii[20]', anatomyDir), # 'passive' tmap
#                sprintf('-c %s/stat_mask_touch_min_rest_ts_concatenated.nii', anatomyDir),
#                sprintf('-expr \u0027step(c)*(a+b)/2\u0027'),
#                sprintf('-prefix %s/stat_active_min_rest_and_passive_min_rest_ts_concatenated.nii', anatomyDir),
#                '-datum FLOAT' )
# system( instr )
# 
# instr <- paste('3dcalc',
#                sprintf('-a %s/stats_concatenated_REML_ts_masked_blurred.nii[17]', anatomyDir),
#                sprintf('-b %s/stats_concatenated_REML_ts_masked_blurred.nii[20]', anatomyDir),
#                sprintf('-c %s/stat_mask_touch_min_rest_REML_ts_concatenated.nii', anatomyDir),
#                sprintf('-expr \u0027step(c)*(a+b)/2\u0027'),
#                sprintf('-prefix %s/stat_active_min_rest_and_passive_min_rest_REML_ts_concatenated.nii', anatomyDir),
#                '-datum FLOAT' )
# system( instr )
# 
# #### ( active + passive ) - ( touch + instructions ) ####
# 
# instr <- paste('3dcalc',
#                sprintf('-a %s/stats_concatenated_ts_masked_blurred.nii[17]', anatomyDir), # 'active' tmap
#                sprintf('-b %s/stats_concatenated_ts_masked_blurred.nii[20]', anatomyDir), # 'passive' tmap
#                sprintf('-c %s/stats_concatenated_ts_masked_blurred.nii[23]', anatomyDir), # 'touch' tmap
#                sprintf('-d %s/stats_concatenated_ts_masked_blurred.nii[26]', anatomyDir), # 'instructions' tmap
#                sprintf('-expr \u0027 (a+b)/2 - (b+c)/2 \u0027'),
#                sprintf('-prefix %s/stat_active_plus_passive_min_touch_plus_instructions_ts_concatenated.nii', anatomyDir),
#                '-datum FLOAT' )
# system( instr )
# 
# instr <- paste('3dcalc',
#                sprintf('-a %s/stats_concatenated_REML_ts_masked_blurred.nii[17]', anatomyDir), # 'active' tmap
#                sprintf('-b %s/stats_concatenated_REML_ts_masked_blurred.nii[20]', anatomyDir), # 'passive' tmap
#                sprintf('-c %s/stats_concatenated_REML_ts_masked_blurred.nii[23]', anatomyDir), # 'touch' tmap
#                sprintf('-d %s/stats_concatenated_REML_ts_masked_blurred.nii[26]', anatomyDir), # 'instructions' tmap
#                sprintf('-expr \u0027 (a+b)/2 - (b+c)/2 \u0027'),
#                sprintf('-prefix %s/stat_active_plus_passive_min_touch_plus_instructions_REML_ts_concatenated.nii', anatomyDir),
#                '-datum FLOAT' )
# system( instr )
# 
# #### ( active - passive ) ####
# 
# instr <- paste('3dcalc',
#                sprintf('-a %s/stats_concatenated_ts_masked_blurred.nii[17]', anatomyDir), # 'active' tmap
#                sprintf('-b %s/stats_concatenated_ts_masked_blurred.nii[20]', anatomyDir), # 'passive' tmap
#                sprintf('-expr \u0027 (a-b) \u0027'),
#                sprintf('-prefix %s/stat_active_min_passive_ts_concatenated.nii', anatomyDir),
#                '-datum FLOAT' )
# system( instr )
# 
# instr <- paste('3dcalc',
#                sprintf('-a %s/stats_concatenated_REML_ts_masked_blurred.nii[17]', anatomyDir), # 'active' tmap
#                sprintf('-b %s/stats_concatenated_REML_ts_masked_blurred.nii[20]', anatomyDir), # 'passive' tmap
#                sprintf('-expr \u0027 (a-b) \u0027'),
#                sprintf('-prefix %s/stat_active_min_passive_REML_ts_concatenated.nii', anatomyDir),
#                '-datum FLOAT' )
# system( instr )
# 


# #### glm on average ts approach ####
# if ( dir.exists('test.results/') ) { system('rm -R test.results/') }
# if ( file.exists('proc.test') ) { system('rm proc.test') }
# if ( file.exists('output.proc.test') ) { system('rm output.proc.test') }
# rString01 <- paste( replicate(4, "\u0027BLOCK(25,1)\u0027"), collapse = " " )
# rString02 <- "\u0027BLOCK(5,1)\u0027"
# 
# instr <- paste('afni_proc.py',
#                '-subj_id test',
#                '-dsets meanTsFolder/meanTs.nii',
#                '-blocks scale regress',
#                '-regress_stim_times stims/*.1D',
#                '-regress_local_times',
#                '-regress_polort 3',
#                '-regress_stim_labels active passive touch rest instructions',
#                sprintf('-regress_opts_3dD'),
#                sprintf('-gltsym \u0027SYM: +active -rest\u0027'),
#                sprintf('-glt_label 1 active-rest'),
#                sprintf('-gltsym \u0027SYM: +passive -rest\u0027'),
#                sprintf('-glt_label 2 passive-rest'),
#                sprintf('-gltsym \u0027SYM: +touch -rest\u0027'),
#                sprintf('-glt_label 3 touch-rest'),
#                sprintf('-gltsym \u0027SYM: +instructions -rest\u0027'),
#                sprintf('-glt_label 4 instructions-rest'),
#                sprintf('-regress_basis_multi %s %s', rString01, rString02 ),
#                '-regress_opts_3dD -jobs 8',
#                '-regress_compute_fitts',
#                '-regress_reml_exec',
#                '-regress_run_clustsim no',
#                '-execute')
# system( instr ) 
# 
# #### bring glm stats from average ts on anatomy ####
# if ( file.exists( sprintf('%s/anatCopy.nii.gz', anatomyDir ) ) ) { system( sprintf('rm %s/anatCopy.nii.gz', anatomyDir ) ) }
# if ( file.exists( sprintf('%s/stats_average_ts.nii', anatomyDir ) ) ) { system( sprintf('rm %s/stats_average_ts.nii', anatomyDir ) ) }
# if ( file.exists( sprintf('%s/stats_average_ts_masked.nii', anatomyDir ) ) ) { system( sprintf('rm %s/stats_average_ts_masked.nii', anatomyDir ) ) }
# if ( file.exists( sprintf('%s/stats_average_ts_masked_blurred.nii', anatomyDir ) ) ) { system( sprintf('rm %s/stats_average_ts_masked_blurred.nii', anatomyDir ) ) }
# if ( file.exists( sprintf('%s/stat_mask_touch_min_rest_ts_average.nii', anatomyDir ) ) ) { system( sprintf('rm %s/stat_mask_touch_min_rest_ts_average.nii', anatomyDir ) ) }
# if ( file.exists( sprintf('%s/stat_active_min_rest_and_passive_min_rest_ts_average.nii', anatomyDir ) ) ) { system( sprintf('rm %s/stat_active_min_rest_and_passive_min_rest_ts_average.nii', anatomyDir ) ) }
# if ( file.exists( sprintf('%s/stat_active_min_rest_and_passive_min_rest_ts_average_blurred.nii', anatomyDir ) ) ) { system( sprintf('rm %s/stat_active_min_rest_and_passive_min_rest_ts_average_blurred.nii', anatomyDir ) ) }
# 
# instr <- sprintf('cp coregistration/anatCopy.nii.gz %s/anatCopy.nii.gz', anatomyDir ); 
# system( instr )
# instr <- sprintf('cp coregistration/singleShot.nii.gz %s/singleShot.nii.gz', anatomyDir ); 
# system( instr )
# instr <- paste( '3dAllineate',
#                 sprintf('-prefix %s/stats_average_ts.nii', anatomyDir),
#                 '-1Dmatrix_apply coregistration/coregMat.1D',
#                 '-final NN',
#                 '-input test.results/stats.test+orig',
#                 sprintf( '-master %s/anatCopy.nii.gz', anatomyDir ) )
# system( instr )
# instr <- paste('3dcalc',
#                sprintf('-a %s/anatCopy.nii.gz', anatomyDir),
#                sprintf('-b %s/stats_average_ts.nii', anatomyDir),
#                sprintf('-expr \u0027step(a)*b\u0027'),
#                sprintf('-prefix %s/stats_average_ts_masked.nii', anatomyDir),
#                '-datum FLOAT' )
# system( instr )
# 
# instr <- paste('3dmerge', 
#                '-1blur_fwhm 1.0',
#                '-doall',
#                sprintf('-prefix %s/stats_average_ts_masked_blurred.nii', anatomyDir),
#                sprintf('%s/stats_average_ts_masked.nii', anatomyDir) )
# system( instr )
# 
# # conjunction analysis from average ts #
# instr <- paste('3dcalc', # mask keeping only those voxels NOT responding to touch-rest,nor to instructions-rest, within the anatomy
#                sprintf('-a %s/anatCopy.nii.gz', anatomyDir),
#                sprintf('-b %s/stats_average_ts_masked_blurred.nii[23]', anatomyDir),
#                sprintf('-c %s/stats_average_ts_masked_blurred.nii[26]', anatomyDir),
#                sprintf('-expr \u0027step(a)*not( step(abs(b)-4) )*not( step(abs(c)-4) )\u0027'),
#                sprintf('-prefix %s/stat_mask_touch_min_rest_ts_average.nii', anatomyDir),
#                '-datum FLOAT' )
# system( instr )
# 
# instr <- paste('3dcalc',
#                sprintf('-a %s/stats_average_ts_masked_blurred.nii[17]', anatomyDir),
#                sprintf('-b %s/stats_average_ts_masked_blurred.nii[20]', anatomyDir),
#                sprintf('-c %s/stat_mask_touch_min_rest_ts_average.nii', anatomyDir),
#                sprintf('-expr \u0027step(c)*(a+b)/2\u0027'),
#                sprintf('-prefix %s/stat_active_min_rest_and_passive_min_rest_ts_average.nii', anatomyDir),
#                '-datum FLOAT' )
# system( instr )

# #### slice time correction ####
# instr <- 'timeSliceCorrection.sh EPI/ *.nii EPI_jsons/ EPI_timeSlice/'; system( instr )
# 
# 
# #### tsnr, standard shimming ####
# system('rm tsnr_standard_shim.nii.gz')
# instr <- '3dTstat -tsnr -prefix tsnr_standard_shim.nii.gz motionCorrect_topUp_Epi/pb.GN19NE455_DanielaTest_MB_GRE_EPI_15VOL_Brainstem_50_slices_test_20220715114301_14.volreg+orig'; system( instr )
# 
# #### tsnr, spm shimming ####
# system('rm tsnr_spm_shim.nii.gz')
# instr <- '3dTstat -tsnr -prefix tsnr_spm_shim.nii.gz motionCorrect_topUp_Epi/pb.GN19NE455_DanielaTest_MB_GRE_EPI_15VOL_Brainstem_50_slices_test_20220715114301_35.volreg+orig'; system( instr )
# 
# #### bring tsnr maps into anatomy, standard shim ####
# instr <- sprintf('rm %s/tsnr_standard_shim_anatSpace.nii.gz', anatomyDir)
# system( instr )
# instr <- paste('3dAllineate', 
#                sprintf('-prefix %s/tsnr_standard_shim_anatSpace.nii.gz',anatomyDir),
#                '-1Dmatrix_apply coregistration/coregMat.1D',
#                '-final NN', 
#                '-input tsnr_standard_shim.nii.gz',
#                sprintf( '-master %s/anatCopy.nii.gz', anatomyDir) )
# print( instr )
# system( instr )
# 
# setwd( anatomyDir )
# prefix_out_file <- 'standard_shim'
# filenameIn <- 'tsnr_standard_shim_anatSpace.nii.gz'
# system( sprintf( 'rm std.141.%s_segvals_05.lh.niml.dset', prefix_out_file ) )
# instr <- paste( '3dVol2Surf',                 
#                 '-spec std.141.anat_test_lh.spec',
#                 '-surf_A std.141.lh.white.gii',
#                 '-surf_B std.141.lh.pial.gii',
#                 sprintf('-sv %s', filenameIn ),
#                 sprintf('-grid_parent %s', filenameIn ),
#                 '-map_func ave',
#                 '-f_steps 5',
#                 '-f_index nodes',
#                 sprintf( '-out_niml std.141.%s_segvals_05.lh.niml.dset', prefix_out_file ) )
# system( instr )  
# 
# setwd( anatomyDir )
# prefix_out_file <- 'standard_shim'
# filenameIn <- 'tsnr_standard_shim_anatSpace.nii.gz'
# system( sprintf( 'rm std.141.%s_segvals_05.rh.niml.dset', prefix_out_file ) )
# instr <- paste( '3dVol2Surf',                 
#                 '-spec std.141.anat_test_rh.spec',
#                 '-surf_A std.141.rh.white.gii',
#                 '-surf_B std.141.rh.pial.gii',
#                 sprintf('-sv %s', filenameIn ),
#                 sprintf('-grid_parent %s', filenameIn ),
#                 '-map_func ave',
#                 '-f_steps 5',
#                 '-f_index nodes',
#                 sprintf( '-out_niml std.141.%s_segvals_05.rh.niml.dset', prefix_out_file ) )
# system( instr )  
# 
# ##### spm shimming ####
# setwd( mainDir )
# instr <- sprintf('rm %s/tsnr_spm_shim_anatSpace.nii.gz', anatomyDir)
# system( instr )
# instr <- paste('3dAllineate', 
#                sprintf('-prefix %s/tsnr_spm_shim_anatSpace.nii.gz',anatomyDir),
#                '-1Dmatrix_apply coregistration/coregMat.1D',
#                '-final NN', 
#                '-input tsnr_spm_shim.nii.gz',
#                sprintf( '-master %s/anatCopy.nii.gz', anatomyDir) )
# print( instr )
# system( instr )
# 
# setwd( anatomyDir )
# prefix_out_file <- 'spm_shim'
# filenameIn <- 'tsnr_spm_shim_anatSpace.nii.gz'
# system( sprintf( 'rm std.141.%s_segvals_05.lh.niml.dset', prefix_out_file ) )
# instr <- paste( '3dVol2Surf',                 
#                 '-spec std.141.anat_test_lh.spec',
#                 '-surf_A std.141.lh.white.gii',
#                 '-surf_B std.141.lh.pial.gii',
#                 sprintf('-sv %s', filenameIn ),
#                 sprintf('-grid_parent %s', filenameIn ),
#                 '-map_func ave',
#                 '-f_steps 5',
#                 '-f_index nodes',
#                 sprintf( '-out_niml std.141.%s_segvals_05.lh.niml.dset', prefix_out_file ) )
# system( instr )  
# 
# setwd( anatomyDir )
# prefix_out_file <- 'spm_shim'
# filenameIn <- 'tsnr_spm_shim_anatSpace.nii.gz'
# system( sprintf( 'rm std.141.%s_segvals_05.rh.niml.dset', prefix_out_file ) )
# instr <- paste( '3dVol2Surf',                 
#                 '-spec std.141.anat_test_rh.spec',
#                 '-surf_A std.141.rh.white.gii',
#                 '-surf_B std.141.rh.pial.gii',
#                 sprintf('-sv %s', filenameIn ),
#                 sprintf('-grid_parent %s', filenameIn ),
#                 '-map_func ave',
#                 '-f_steps 5',
#                 '-f_index nodes',
#                 sprintf( '-out_niml std.141.%s_segvals_05.rh.niml.dset', prefix_out_file ) )
# system( instr )  
# 
# #### bring epi into anat ####
# setwd( mainDir )
# instr <- sprintf('rm %s/epi_anatSpace.nii.gz', anatomyDir)
# system( instr )
# instr <- paste('3dAllineate', 
#                sprintf('-prefix %s/epi_anatSpace.nii.gz',anatomyDir),
#                '-1Dmatrix_apply coregistration/coregMat.1D',
#                '-final lin', 
#                '-input meanEpi4Coreg.nii.gz',
#                sprintf( '-master %s/anatCopy.nii.gz', anatomyDir) )
# print( instr )
# system( instr )
# 
# #### GLM saccades ####
# setwd( mainDir )
# system('rm -R glm_saccades.results')
# system('rm *glm_saccades*')
# instr <- paste('afni_proc.py',
#   '-subj_id glm_saccades',
#   '-dsets EPI_timeSlice/GN19NE455_DanielaTest_MB_GRE_EPI_15VOL_Brainstem_50_slices_test_20220715114301_35.nii',
#   '-blocks despike volreg mask scale regress',
#   '-regress_stim_times stimTimes/saccadeTimes.1D',
#   '-regress_stim_labels sacc_on',
#   sprintf( '-regress_basis \u0027BLOCK(15,1)\u0027' ) )
# system( instr )
# system('tcsh -xef proc.glm_saccades 2>&1 | tee output.proc.glm_saccades')
# 
# #### GLM finger/eye ####
# setwd( mainDir )
# system('rm -R glm_movement.results')
# system('rm *glm_movement*')
# instr <- paste('afni_proc.py',
#                '-subj_id glm_movement',
#                '-dsets EPI_timeSlice/GN19NE455_DanielaTest_MB_GRE_EPI_15VOL_Brainstem_50_slices_test_20220715114301_38.nii',
#                '-blocks despike volreg mask scale regress',
#                '-regress_stim_times stimTimes/movementTimes.1D',
#                '-regress_stim_labels sacc_on',
#                sprintf( '-regress_basis \u0027BLOCK(15,1)\u0027' ) )
# system( instr )
# system('tcsh -xef proc.glm_movement 2>&1 | tee output.proc.glm_movement')
# 
# 
# #### atlas into anatomy, lh hemi ####
# setwd( anatomyDir )
# prefix_out_file <- 'princeton_rois_lh.nii.gz'
# filenameIn <- 'maxprob_surf_lh.1D.dset'
# system( sprintf( 'rm %s', prefix_out_file ) )
# instr <- paste( '3dSurf2Vol',                 
#                 '-spec std.141.anat_test_lh.spec',
#                 '-surf_A std.141.lh.white.gii',
#                 '-surf_B std.141.lh.pial.gii',
#                 '-sv anatCopy.nii.gz',
#                 '-grid_parent anatCopy.nii.gz',
#                 '-sdata maxprob_surf_lh.1D.dset',
#                 '-map_func mode',
#                 '-f_steps 5',
#                 '-f_index nodes',
#                 sprintf( '-prefix %s', prefix_out_file ) )
# system( instr )  
# 
# 
# #### atlas into anatomy, rh hemi ####
# setwd( anatomyDir )
# prefix_out_file <- 'princeton_rois_rh.nii.gz'
# filenameIn <- 'maxprob_surf_rh.1D.dset'
# system( sprintf( 'rm %s', prefix_out_file ) )
# instr <- paste( '3dSurf2Vol',                 
#                 '-spec std.141.anat_test_rh.spec',
#                 '-surf_A std.141.rh.white.gii',
#                 '-surf_B std.141.rh.pial.gii',
#                 '-sv anatCopy.nii.gz',
#                 '-grid_parent anatCopy.nii.gz',
#                 '-sdata maxprob_surf_rh.1D.dset',
#                 '-map_func mode',
#                 '-f_steps 5',
#                 '-f_index nodes',
#                 sprintf( '-prefix %s', prefix_out_file ) )
# system( instr )  
# 
# #### atlas into anatomy, Glasser, lh hemi ####
# setwd( anatomyDir )
# prefix_out_file <- 'Glasser_rois_lh.nii.gz'
# filenameIn <- 'lh.std.141.Glasser_HCP.niml.dset'
# system( sprintf( 'rm %s', prefix_out_file ) )
# instr <- paste( '3dSurf2Vol',                 
#                 '-spec std.141.anat_test_lh.spec',
#                 '-surf_A std.141.lh.white.gii',
#                 '-surf_B std.141.lh.pial.gii',
#                 '-sv anatCopy.nii.gz',
#                 '-grid_parent anatCopy.nii.gz',
#                 sprintf('-sdata %s',filenameIn),
#                 '-map_func mode',
#                 '-f_steps 5',
#                 '-f_index nodes',
#                 sprintf( '-prefix %s', prefix_out_file ) )
# system( instr )  
# 
# 
# #### atlas into anatomy, Glasser, rh hemi ####
# setwd( anatomyDir )
# prefix_out_file <- 'Glasser_rois_rh.nii.gz'
# filenameIn <- 'rh.std.141.Glasser_HCP.niml.dset'
# system( sprintf( 'rm %s', prefix_out_file ) )
# instr <- paste( '3dSurf2Vol',                 
#                 '-spec std.141.anat_test_rh.spec',
#                 '-surf_A std.141.rh.white.gii',
#                 '-surf_B std.141.rh.pial.gii',
#                 '-sv anatCopy.nii.gz',
#                 '-grid_parent anatCopy.nii.gz',
#                 sprintf('-sdata %s',filenameIn),
#                 '-map_func mode',
#                 '-f_steps 5',
#                 '-f_index nodes',
#                 sprintf( '-prefix %s', prefix_out_file ) )
# system( instr )  
# 
# 
# #### atlas, ROI and tsnr into dataframe ####
# afniDir <- Sys.getenv(x='AFNI_INSTALLDIR')
# afniSurfaces <- Sys.getenv(x='AFNI_TOOLBOXDIRSURFACES')
# source( sprintf( '%s/AFNIio.R', afniDir ) )
# source( sprintf( '%s/coordinateFromLinearIndex.r', afniSurfaces ) )
# getRoiIdx <- function( roiFileName, volumeIdx ) {
#   roiFile <- read.AFNI( roiFileName )
#   roiVolume <- roiFile$brk[,,,volumeIdx]
#   roiIdxTemp <- unique( array( roiVolume ) )
#   roiIdx <- sort( roiIdxTemp[ roiIdxTemp>0 ] )
#   for (i in roiIdx) {
#     idxVolTemp <- which( roiVolume==roiIdx[i] )
#     roiLabel <- rep(i,length(idxVolTemp))
#     if (i==1) { dfOut <- data.frame( roiLabel, idxVolTemp ) }
#     if (i>1) {
#       dfOutTemp <- data.frame( roiLabel, idxVolTemp )
#       dfOut <- rbind( dfOut, dfOutTemp )
#     }
#   }
#   return( dfOut )
# }
# getRoiData <- function( inputFilename, roisDF, volumeIdx ) {
#   dataFile <- read.AFNI( inputFilename )
#   dataVolume <- dataFile$brk[,,,volumeIdx]
#   roiIdx <- unique( roisDF$roiLabel )
#   roisDF$volumeValue <- rep(9999, dim(roisDF)[1])
#   for (i in roiIdx) {
#     selectedRoiIdx <- roisDF$idxVolTemp[ roisDF$roiLabel==i ]
#     roisDF$volumeValue[ roisDF$roiLabel==i ] <-  dataVolume[ selectedRoiIdx ]
#   }
#   return( roisDF )
# }
# 
# setwd(anatomyDir)
# dfRoisLh <- getRoiIdx( 'princeton_rois_lh.nii.gz', 2 )
# dfRoisRh <- getRoiIdx( 'princeton_rois_rh.nii.gz', 2 )
# dfRoisBrainstem <- getRoiIdx( 'subCorticalRois+orig', 1 )
# dfRoisLh_Glasser <- getRoiIdx( 'Glasser_rois_lh.nii.gz', 1 )
# dfRoisRh_Glasser <- getRoiIdx( 'Glasser_rois_rh.nii.gz', 1 )
# 
# dfRois_tsnr_standard_lh <- getRoiData( 'tsnr_standard_shim_anatSpace.nii.gz', dfRoisLh, 1  )
# dfRois_tsnr_standard_rh <- getRoiData( 'tsnr_standard_shim_anatSpace.nii.gz', dfRoisRh, 1  )
# dfRois_tsnr_standard_brainstem <- getRoiData( 'tsnr_standard_shim_anatSpace.nii.gz', dfRoisBrainstem, 1  )
# dfRois_tsnr_Glasser_lh_standard <- getRoiData( 'tsnr_standard_shim_anatSpace.nii.gz', dfRoisLh_Glasser, 1  )
# dfRois_tsnr_Glasser_rh_standard <- getRoiData( 'tsnr_standard_shim_anatSpace.nii.gz', dfRoisRh_Glasser, 1  )
# 
# dfRois_tsnr_spm_lh <- getRoiData( 'tsnr_spm_shim_anatSpace.nii.gz', dfRoisLh, 1  )
# dfRois_tsnr_spm_rh <- getRoiData( 'tsnr_spm_shim_anatSpace.nii.gz', dfRoisRh, 1  )
# dfRois_tsnr_spm_brainstem <- getRoiData( 'tsnr_spm_shim_anatSpace.nii.gz', dfRoisBrainstem, 1  )
# dfRois_tsnr_Glasser_lh_spm <- getRoiData( 'tsnr_spm_shim_anatSpace.nii.gz', dfRoisLh_Glasser, 1  )
# dfRois_tsnr_Glasser_rh_spm <- getRoiData( 'tsnr_spm_shim_anatSpace.nii.gz', dfRoisRh_Glasser, 1  )
# 
# lhStandardDF <- aggregate( volumeValue ~ roiLabel, median, data=dfRois_tsnr_standard_lh )
# rhStandardDF <- aggregate( volumeValue ~ roiLabel, median, data=dfRois_tsnr_standard_rh )
# brainstemStandardDF <- aggregate( volumeValue ~ roiLabel, median, data=dfRois_tsnr_standard_brainstem )
# lhStandardDF_sd <- aggregate( volumeValue ~ roiLabel, sd, data=dfRois_tsnr_standard_lh )
# rhStandardDF_sd <- aggregate( volumeValue ~ roiLabel, sd, data=dfRois_tsnr_standard_rh )
# brainstemStandardDF_sd <- aggregate( volumeValue ~ roiLabel, sd, data=dfRois_tsnr_standard_brainstem )
# lhStandardDF_Glasser <- aggregate( volumeValue ~ roiLabel, median, data=dfRois_tsnr_Glasser_lh )
# rhStandardDF_Glasser <- aggregate( volumeValue ~ roiLabel, median, data=dfRois_tsnr_Glasser_rh )
# lhStandardDF_Glasser_sd <- aggregate( volumeValue ~ roiLabel, sd, data=dfRois_tsnr_Glasser_lh )
# rhStandardDF_Glasser_sd <- aggregate( volumeValue ~ roiLabel, sd, data=dfRois_tsnr_Glasser_rh )
# 
# lhSpmDF <- aggregate( volumeValue ~ roiLabel, median, data=dfRois_tsnr_spm_lh )
# rhSpmDF <- aggregate( volumeValue ~ roiLabel, median, data=dfRois_tsnr_spm_rh )
# brainstemSpmDF <- aggregate( volumeValue ~ roiLabel, median, data=dfRois_tsnr_spm_brainstem )
# lhSpmDF_sd <- aggregate( volumeValue ~ roiLabel, sd, data=dfRois_tsnr_spm_lh )
# rhSpmDF_sd <- aggregate( volumeValue ~ roiLabel, sd, data=dfRois_tsnr_spm_rh )
# brainstemSpmDF_sd <- aggregate( volumeValue ~ roiLabel, sd, data=dfRois_tsnr_spm_brainstem )
# lhStandardDF_Glasser_spm <- aggregate( volumeValue ~ roiLabel, median, data=dfRois_tsnr_spm_brainstem )
# rhStandardDF_Glasser_spm <- aggregate( volumeValue ~ roiLabel, median, data=dfRois_tsnr_spm_brainstem )
# lhStandardDF_Glasser_spm_sd <- aggregate( volumeValue ~ roiLabel, sd, data=dfRois_tsnr_spm_brainstem )
# rhStandardDF_Glasser_spm_sd <- aggregate( volumeValue ~ roiLabel, sd, data=dfRois_tsnr_spm_brainstem )
# 
# #### plots standard shimming ####
# x11( width=6, height=6)
# plotDf <- subset( dfRois_tsnr_standard_lh, roiLabel==1 | roiLabel==2 | roiLabel==3 )
# boxplot( plotDf$volumeValue ~ plotDf$roiLabel, las=1, frame=FALSE, outline=FALSE, boxwex=0.25,
#          at=c(1-0.15,2-0.15,3-0.15), xlim=c(0.5,3.5), axes=FALSE, col='lightblue', ylim=c(0,150),
#          xlab='ROI', ylab='tSNR', cex.lab=2, lwd=2 )
# plotDf <- subset( dfRois_tsnr_standard_rh, roiLabel==1 | roiLabel==2 | roiLabel==3 )
# boxplot( plotDf$volumeValue ~ plotDf$roiLabel, las=1, frame=FALSE, outline=FALSE, boxwex=0.25,
#          at=c(1+0.15,2+0.15,3+0.15), xlim=c(0,4), axes=FALSE, col='orange', add=TRUE, lwd=2 )
# axis(1, seq(1,3), c('V1','V2','V3'), cex.axis=2 )
# axis(2, seq(0,200,50), cex.axis=2, las=1 )
# dev.copy2pdf( file=sprintf('%s/standardShim_earlyVC.pdf', figuresDir) )
# graphics.off()
# 
# x11( width=6, height=6)
# plotDf <- subset( dfRois_tsnr_standard_brainstem, roiLabel==1 | roiLabel==2 | roiLabel==3 )
# boxplot( plotDf$volumeValue ~ plotDf$roiLabel, las=1, frame=FALSE, outline=FALSE, boxwex=0.75,
#          at=c(1,2,3), xlim=c(0.5,3.5), axes=FALSE, col=c('red','lightblue','orange'), ylim=c(0,150),
#          xlab='ROI', ylab='tSNR', cex.lab=2, lwd=2 )
# axis(1, seq(1,3), c('V1','V2','V3'), cex.axis=2 )
# axis(2, seq(0,200,50), cex.axis=2, las=1 )
# dev.copy2pdf( file=sprintf('%s/standardShim_brainstem.pdf', figuresDir) )
# graphics.off()
# 
# #### plots spm shimming ####
# x11( width=6, height=6)
# plotDf <- subset( dfRois_tsnr_spm_lh, roiLabel==1 | roiLabel==2 | roiLabel==3 )
# boxplot( plotDf$volumeValue ~ plotDf$roiLabel, las=1, frame=FALSE, outline=FALSE, boxwex=0.25,
#          at=c(1-0.15,2-0.15,3-0.15), xlim=c(0.5,3.5), axes=FALSE, col='lightblue', ylim=c(0,150),
#          xlab='ROI', ylab='tSNR', cex.lab=2, lwd=2 )
# plotDf <- subset( dfRois_tsnr_spm_rh, roiLabel==1 | roiLabel==2 | roiLabel==3 )
# boxplot( plotDf$volumeValue ~ plotDf$roiLabel, las=1, frame=FALSE, outline=FALSE, boxwex=0.25,
#          at=c(1+0.15,2+0.15,3+0.15), xlim=c(0,4), axes=FALSE, col='orange', add=TRUE, lwd=2 )
# axis(1, seq(1,3), c('V1','V2','V3'), cex.axis=2 )
# axis(2, seq(0,200,50), cex.axis=2, las=1 )
# dev.copy2pdf( file=sprintf('%s/spmShim_earlyVC.pdf', figuresDir) )
# graphics.off()
# 
# x11( width=6, height=6)
# plotDf <- subset( dfRois_tsnr_spm_brainstem, roiLabel==1 | roiLabel==2 | roiLabel==3 )
# boxplot( plotDf$volumeValue ~ plotDf$roiLabel, las=1, frame=FALSE, outline=FALSE, boxwex=0.75,
#          at=c(1,2,3), xlim=c(0.5,3.5), axes=FALSE, col=c('red','lightblue','orange'), ylim=c(0,150),
#          xlab='ROI', ylab='tSNR', cex.lab=2, lwd=2 )
# axis(1, seq(1,3), c('V1','V2','V3'), cex.axis=2 )
# axis(2, seq(0,200,50), cex.axis=2, las=1 )
# dev.copy2pdf( file=sprintf('%s/spmShim_brainstem.pdf', figuresDir) )
# graphics.off()
# 
# #### plots standard shimming, evc and hippocampus ####
# x11( width=11, height=6)
# plotDf <- subset( dfRois_tsnr_Glasser_lh_standard, roiLabel==1 | roiLabel==2 | roiLabel==3 | roiLabel==120 | roiLabel==119 | roiLabel==126 | roiLabel==155 | roiLabel==127 | roiLabel==118 )
# boxplot( plotDf$volumeValue ~ plotDf$roiLabel, las=1, frame=FALSE, outline=FALSE, boxwex=0.25,
#          at=c(1-0.15,2-0.15,3-0.15,4-0.15,5-0.15,6-0.15,7-0.15,8-0.15,9-0.15), xlim=c(0.5,9.5), axes=FALSE, col='lightblue', ylim=c(0,150),
#          xlab='ROI', ylab='tSNR', cex.lab=2, lwd=2 )
# plotDf <- subset( dfRois_tsnr_Glasser_rh_standard, roiLabel==1 | roiLabel==2 | roiLabel==3 | roiLabel==120 | roiLabel==119 | roiLabel==126 | roiLabel==155 | roiLabel==127 | roiLabel==118 )
# boxplot( plotDf$volumeValue ~ plotDf$roiLabel, las=1, frame=FALSE, outline=FALSE, boxwex=0.25,
#          at=c(1+0.15,2+0.15,3+0.15,4+0.15,5+0.15,6+0.15,7+0.15,8+0.15,9+0.15), xlim=c(0.5,9.5), axes=FALSE, col='orange', add=TRUE, lwd=2 )
# axis(1, seq(1,9), c('V1','V2','V3','H','PRE','PH1','PH2','PH3','ETH'), cex.axis=2 )
# axis(2, seq(0,200,50), cex.axis=2, las=1 )
# dev.copy2pdf( file=sprintf('%s/standardShim_earlyVC_hippocampus.pdf', figuresDir) )
# graphics.off()
# 
# #### plots sp, shimming, evc and hippocampus ####
# x11( width=11, height=6)
# plotDf <- subset( dfRois_tsnr_Glasser_lh_spm, roiLabel==1 | roiLabel==2 | roiLabel==3 | roiLabel==120 | roiLabel==119 | roiLabel==126 | roiLabel==155 | roiLabel==127 | roiLabel==118 )
# boxplot( plotDf$volumeValue ~ plotDf$roiLabel, las=1, frame=FALSE, outline=FALSE, boxwex=0.25,
#          at=c(1-0.15,2-0.15,3-0.15,4-0.15,5-0.15,6-0.15,7-0.15,8-0.15,9-0.15), xlim=c(0.5,9.5), axes=FALSE, col='lightblue', ylim=c(0,150),
#          xlab='ROI', ylab='tSNR', cex.lab=2, lwd=2 )
# plotDf <- subset( dfRois_tsnr_Glasser_rh_spm, roiLabel==1 | roiLabel==2 | roiLabel==3 | roiLabel==120 | roiLabel==119 | roiLabel==126 | roiLabel==155 | roiLabel==127 | roiLabel==118 )
# boxplot( plotDf$volumeValue ~ plotDf$roiLabel, las=1, frame=FALSE, outline=FALSE, boxwex=0.25,
#          at=c(1+0.15,2+0.15,3+0.15,4+0.15,5+0.15,6+0.15,7+0.15,8+0.15,9+0.15), xlim=c(0.5,9.5), axes=FALSE, col='orange', add=TRUE, lwd=2 )
# axis(1, seq(1,9), c('V1','V2','V3','H','PRE','PH1','PH2','PH3','ETH'), cex.axis=2 )
# axis(2, seq(0,200,50), cex.axis=2, las=1 )
# dev.copy2pdf( file=sprintf('%s/spmShim_earlyVC_hippocampus.pdf', figuresDir) )
# graphics.off()
