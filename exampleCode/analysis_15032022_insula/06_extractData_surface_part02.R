rm( list=ls() )
afniDir <- Sys.getenv(x='AFNI_INSTALLDIR')
afniSurfaces <- Sys.getenv(x='AFNI_TOOLBOXDIRSURFACES')
source( sprintf( '%s/AFNIio.R', afniDir ) )
source( sprintf( '%s/coordinateFromLinearIndex.r', afniSurfaces ) )

mainDir <- '/analyse/Project0165/Fracasso_Anatomy/AHEAD_database/10007840'
partDir <- 'part2/'

profileFunction <- function( filenameIn, hemi, prefix_out_file ) {
  
  if (hemi=='LH') {
    # left hemisphere .niml
    system( sprintf( 'rm std.141.%s_segvals_10.lh.niml.dset', prefix_out_file ) )
    instr <- paste( '3dVol2Surf',                 
                    '-spec std.141.AHEAD_test_lh.spec',
                    '-surf_A std.141.lh.white.gii',
                    '-surf_B std.141.lh.pial.gii',
                    sprintf('-sv %s', filenameIn ),
                    sprintf('-grid_parent %s', filenameIn ),
                    '-map_func seg_vals',
                    '-f_steps 10',
                    '-f_index nodes',
                    sprintf( '-out_niml std.141.%s_segvals_10.lh.niml.dset', prefix_out_file ) )
    system( instr )  
    
    # left hemisphere .1D
    system( sprintf( 'rm std.141.%s_segvals_10.lh.1D.dset', prefix_out_file ) )
    instr <- paste( '3dVol2Surf',                 
                    '-spec std.141.AHEAD_test_lh.spec',
                    '-surf_A std.141.lh.white.gii',
                    '-surf_B std.141.lh.pial.gii',
                    sprintf('-sv %s', filenameIn ),
                    sprintf('-grid_parent %s', filenameIn ),
                    '-map_func seg_vals',
                    '-f_steps 10',
                    '-f_index nodes',
                    sprintf( '-out_1D std.141.%s_segvals_10.lh.1D.dset', prefix_out_file ) )
    system( instr )  
  }
  
  if (hemi=='RH') {
    # right hemisphere .niml
    system( sprintf( 'rm std.141.%s_segvals_10.rh.niml.dset', prefix_out_file ) )
    instr <- paste( '3dVol2Surf',                 
                    '-spec std.141.AHEAD_test_rh.spec',
                    '-surf_A std.141.rh.white.gii',
                    '-surf_B std.141.rh.pial.gii',
                    sprintf('-sv %s', filenameIn ),
                    sprintf('-grid_parent %s', filenameIn ),
                    '-map_func seg_vals',
                    '-f_steps 10',
                    '-f_index nodes',
                    sprintf( '-out_niml std.141.%s_segvals_10.rh.niml.dset', prefix_out_file ) )
    system( instr )  
    
    # right hemisphere .1D
    system( sprintf( 'rm std.141.%s_segvals_10.rh.1D.dset', prefix_out_file ) )
    instr <- paste( '3dVol2Surf',                 
                    '-spec std.141.AHEAD_test_rh.spec',
                    '-surf_A std.141.rh.white.gii',
                    '-surf_B std.141.rh.pial.gii',
                    sprintf('-sv %s', filenameIn ),
                    sprintf('-grid_parent %s', filenameIn ),
                    '-map_func seg_vals',
                    '-f_steps 10',
                    '-f_index nodes',
                    sprintf( '-out_1D std.141.%s_segvals_10.rh.1D.dset', prefix_out_file ) )
    system( instr )  
  }
  
}

setwd( mainDir )
setwd( partDir )
getwd()
participantsFolders <- dir()
for ( participantsFoldersIndex in seq( 1, ( length( participantsFolders ) ) ) ) { #seq( 1,length( participantsFolders )
  
  setwd( mainDir )
  setwd( partDir )
  print( '...' )
  print( '...' )
  print( sprintf( '%s...', participantsFolders[ participantsFoldersIndex ] ) )
  print( '...' )
  print( '...' )
  setwd( sprintf( '%s', participantsFolders[ participantsFoldersIndex ] ) )
  setwd( dir()[1] )
  setwd( dir()[1] )
  print( dir() )
  setwd( 'AHEAD_test' )
  setwd( 'SUMA' )
  print( dir() )
  
  # left hemisphere atlas, .1D
  system('rm std.141.aparc+aseg_REN_all_midpoint.lh.1D.dset')
  instr <- paste( '3dVol2Surf',                 
                  '-spec std.141.AHEAD_test_lh.spec',
                  '-surf_A std.141.lh.white.gii',
                  '-surf_B std.141.lh.pial.gii',
                  '-sv aparc+aseg_REN_all.nii.gz',
                  '-grid_parent aparc+aseg_REN_all.nii.gz',
                  '-map_func midpoint',
                  '-f_index nodes',
                  '-out_1D std.141.aparc+aseg_REN_all_midpoint.lh.1D.dset' )
  system( instr )  
  
  # right hemisphere atlas, .1D
  system('rm std.141.aparc+aseg_REN_all_midpoint.rh.1D.dset')
  instr <- paste( '3dVol2Surf',                 
                  '-spec std.141.AHEAD_test_rh.spec',
                  '-surf_A std.141.rh.white.gii',
                  '-surf_B std.141.rh.pial.gii',
                  '-sv aparc+aseg_REN_all.nii.gz',
                  '-grid_parent aparc+aseg_REN_all.nii.gz',
                  '-map_func midpoint',
                  '-f_index nodes',
                  '-out_1D std.141.aparc+aseg_REN_all_midpoint.rh.1D.dset' )
  system( instr )  
  
  # curvature, 1D
  system('rm std.141.lh.curv_ppp_.1D.dset std.141.rh.curv_ppp_.1D.dset')
  instr <- 'ConvertDset -o_1D -input std.141.lh.curv.niml.dset -prepend_node_index_1D -prefix std.141.lh.curv_ppp_.1D.dset'; system( instr )
  instr <- 'ConvertDset -o_1D -input std.141.rh.curv.niml.dset -prepend_node_index_1D -prefix std.141.rh.curv_ppp_.1D.dset'; system( instr )
  
  # thickness, 1D
  system('rm std.141.lh.thickness_ppp_.1D.dset std.141.rh.thickness_ppp_.1D.dset')
  instr <- 'ConvertDset -o_1D -input std.141.lh.thickness.niml.dset -prepend_node_index_1D -prefix std.141.lh.thickness_ppp_.1D.dset'; system( instr )
  instr <- 'ConvertDset -o_1D -input std.141.rh.thickness.niml.dset -prepend_node_index_1D -prefix std.141.rh.thickness_ppp_.1D.dset'; system( instr )

  # sulc, 1D
  system('rm std.141.lh.sulc_ppp_.1D.dset std.141.rh.sulc_ppp_.1D.dset')
  instr <- 'ConvertDset -o_1D -input std.141.lh.sulc.niml.dset -prepend_node_index_1D -prefix std.141.lh.sulc_ppp_.1D.dset'; system( instr )
  instr <- 'ConvertDset -o_1D -input std.141.rh.sulc.niml.dset -prepend_node_index_1D -prefix std.141.rh.sulc_ppp_.1D.dset'; system( instr )
  
  subPrefix <- participantsFolders[ participantsFoldersIndex ]
  
  profileFunction( sprintf( '%s_ses-1_acq-wb_mod-t1w_orient-std_brain.nii.gz', subPrefix ), 'LH', 't1w' )
  profileFunction( sprintf( '%s_ses-1_acq-wb_mod-t1w_orient-std_brain.nii.gz', subPrefix ), 'RH', 't1w' )
  
  profileFunction( sprintf( '%s_ses-1_acq-wb_mod-t1map_orient-std_brain.nii.gz', subPrefix ), 'LH', 't1map' )
  profileFunction( sprintf( '%s_ses-1_acq-wb_mod-t1map_orient-std_brain.nii.gz', subPrefix ), 'RH', 't1map' )
  
  profileFunction( sprintf( '%s_ses-1_acq-wb_mod-r1map_orient-std_brain.nii.gz', subPrefix ), 'LH', 'r1map' )
  profileFunction( sprintf( '%s_ses-1_acq-wb_mod-r1map_orient-std_brain.nii.gz', subPrefix ), 'RH', 'r1map' )
  
  profileFunction( sprintf( '%s_ses-1_acq-wb_mod-qsm_orient-std_brain.nii.gz', subPrefix ), 'LH', 'qsm' )
  profileFunction( sprintf( '%s_ses-1_acq-wb_mod-qsm_orient-std_brain.nii.gz', subPrefix ), 'RH', 'qsm' )

  profileFunction( sprintf( '%s_ses-1_acq-wb_mod-t2starmap_orient-std_brain.nii.gz', subPrefix ), 'LH', 't2starmap' )
  profileFunction( sprintf( '%s_ses-1_acq-wb_mod-t2starmap_orient-std_brain.nii.gz', subPrefix ), 'RH', 't2starmap' )
  
  profileFunction( sprintf( '%s_ses-1_acq-wb_mod-r2starmap_orient-std_brain.nii.gz', subPrefix ), 'LH', 'r2starmap' )
  profileFunction( sprintf( '%s_ses-1_acq-wb_mod-r2starmap_orient-std_brain.nii.gz', subPrefix ), 'RH', 'r2starmap' )
  
}