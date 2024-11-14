rm(list=ls())
atlasFolder <- '/analyse/Project0226/tests/Daniela_pilotData/suma_MNI152_2009_princetonAtlas'
  
for (nSubj in 1:8) { #nSubj <- 6
  
  if ( nSubj==1 ) {
    mainDir <- '/analyse/Project0226/tests/Daniela_pilotData/DanielaProject03082023_Daniela'
    anatomyDir <- '/analyse/Project0226/tests/Daniela_pilotData/DanielaProject03082023_Daniela/ANATOMY/anat_test/SUMA'
  }
  
  if ( nSubj==2 ) {
    mainDir <- '/analyse/Project0226/tests/Daniela_pilotData/DanielaProject22082023_Holly'
    anatomyDir <- '/analyse/Project0226/tests/Daniela_pilotData/DanielaProject22082023_Holly/ANATOMY/anat_test/SUMA'
  }
    
  if ( nSubj==3 ) {
    mainDir <- '/analyse/Project0226/tests/Daniela_pilotData/DanielaProject22082023_Belinda'
    anatomyDir <- '/analyse/Project0226/tests/Daniela_pilotData/DanielaProject22082023_Belinda/ANATOMY/anat_test/SUMA'
  }
  
  if ( nSubj==4 ) {
    mainDir <- '/analyse/Project0226/tests/Daniela_pilotData/DanielaProject28092023_MTR13'
    anatomyDir <- '/analyse/Project0226/tests/Daniela_pilotData/DanielaProject28092023_MTR13/ANATOMY/anat_test/SUMA'
  }
  
  if ( nSubj==5 ) {
    mainDir <- '/analyse/Project0226/tests/Daniela_pilotData/DanielaProject28092023_NWS30'
    anatomyDir <- '/analyse/Project0226/tests/Daniela_pilotData/DanielaProject28092023_NWS30/ANATOMY/anat_test/SUMA'
  }

  if ( nSubj==6 ) {
    mainDir <- '/analyse/Project0226/tests/Daniela_pilotData/DanielaProject10102023_Rosie'
    anatomyDir <- '/analyse/Project0226/tests/Daniela_pilotData/DanielaProject10102023_Rosie/ANATOMY/anat_test/SUMA'
  }
  
  if ( nSubj==7 ) {
    mainDir <- '/analyse/Project0226/tests/Daniela_pilotData/DanielaProject31102023_LucasEdward'
    anatomyDir <- '/analyse/Project0226/tests/Daniela_pilotData/DanielaProject31102023_LucasEdward/ANATOMY/anat_test/SUMA'
  }

  if ( nSubj==8 ) {
    mainDir <- '/analyse/Project0226/tests/Daniela_pilotData/DanielaProject31102023_Theo'
    anatomyDir <- '/analyse/Project0226/tests/Daniela_pilotData/DanielaProject31102023_Theo/ANATOMY/anat_test/SUMA'
  }

  setwd( mainDir )
  getwd()
  
  subjFilename <- strsplit( mainDir, '[_]')[[1]][3]
  
  setwd( atlasFolder )
  getwd()
  
  # clean up
  if ( file.exists( sprintf( '%s/stats_concatenated_REML_ts_masked_%s.nii', atlasFolder, subjFilename ) ) ) {
    instr <- sprintf( 'rm %s/stats_concatenated_REML_ts_masked_%s.nii', atlasFolder, subjFilename ); system( instr )
  } 
  if ( file.exists( sprintf( '%s/stats_concatenated_ts_masked_%s.nii', atlasFolder, subjFilename ) ) ) {
    instr <- sprintf( 'rm %s/stats_concatenated_ts_masked_%s.nii', atlasFolder, subjFilename ); system( instr )
  } 
  if ( file.exists( sprintf( '%s/stats_concatenated_REML_ts_masked_linear_%s.nii', atlasFolder, subjFilename ) ) ) {
    instr <- sprintf( 'rm %s/stats_concatenated_REML_ts_masked_linear_%s.nii', atlasFolder, subjFilename ); system( instr )
  } 
  if ( file.exists( sprintf( '%s/stats_concatenated_ts_masked_linear_%s.nii', atlasFolder, subjFilename ) ) ) {
    instr <- sprintf( 'rm %s/stats_concatenated_ts_masked_linear_%s.nii', atlasFolder, subjFilename ); system( instr )
  } 
  
  # bring stats REML into anatomy, non-linear
  instr <- paste('3dNwarpApply',
                 sprintf( '-master brain.nii' ),
                 sprintf( '-source %s/stats_concatenated_REML_ts_masked.nii', anatomyDir ),
                 sprintf('-nwarp Qwarp_WARP_%s+tlrc MPRAGE_%s_Xat.1D', subjFilename, subjFilename ),
                 '-dxyz 2 2 2',
                 '-interp NN',
                 '-short')
  print( instr )
  system( instr )
  system( sprintf('mv stats_concatenated_REML_ts_masked_Nwarp.nii stats_concatenated_REML_ts_masked_%s.nii', subjFilename ) )

  # bring stats REML into anatomy, linear
  instr <- paste( '3dAllineate',
                  sprintf( '-prefix stats_concatenated_REML_ts_masked_linear_%s.nii', subjFilename ),
                  sprintf( '-1Dmatrix_apply MPRAGE_%s_Xaff12.1D', subjFilename ), 
                  '-final NN', 
                  '-mast_dxyz 2',
                  sprintf('-input %s/stats_concatenated_REML_ts_masked.nii', anatomyDir ),
                  '-master brain.nii' )
  print( instr )
  system( instr )
  
  # bring stats into anatomy, non-linear
  instr <- paste('3dNwarpApply',
                 sprintf( '-master brain.nii' ),
                 sprintf( '-source %s/stats_concatenated_ts_masked.nii', anatomyDir ),
                 sprintf('-nwarp Qwarp_WARP_%s+tlrc MPRAGE_%s_Xat.1D ', subjFilename, subjFilename ),
                 '-dxyz 2 2 2',
                 '-interp NN',
                 '-short')
  print( instr )
  system( instr )
  system( sprintf('mv stats_concatenated_ts_masked_Nwarp.nii stats_concatenated_ts_masked_%s.nii', subjFilename ) )
  
  # bring stats into anatomy, linear
  instr <- paste( '3dAllineate',
                  sprintf( '-prefix stats_concatenated_ts_masked_linear_%s.nii', subjFilename ),
                  sprintf( '-1Dmatrix_apply MPRAGE_%s_Xaff12.1D', subjFilename ), 
                  '-final NN', 
                  '-mast_dxyz 2',
                  sprintf('-input %s/stats_concatenated_ts_masked.nii', anatomyDir ),
                  '-master brain.nii' )
  print( instr )
  system( instr )
  
  
}

