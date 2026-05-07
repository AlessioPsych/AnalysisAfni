rm( list=ls() )
afniDir <- Sys.getenv(x='AFNI_INSTALLDIR')
afniSurfaces <- Sys.getenv(x='AFNI_TOOLBOXDIRSURFACES')
source( sprintf( '%s/AFNIio.R', afniDir ) )
source( sprintf( '%s/coordinateFromLinearIndex.r', afniSurfaces ) )

saveRoi <- '/analyse/Project0165/Fracasso_Anatomy/AHEAD_database/10007840/results/insula'
appendFilename <- 'LH'

mainDir <- '/analyse/Project0165/Fracasso_Anatomy/AHEAD_database/10007840'
partDir <- 'part1/'

setwd( mainDir )
setwd( partDir )
getwd()
participantsFolders <- dir()
for ( participantsFoldersIndex in 2:4 ) { #seq( 1,length( participantsFolders )
  
  #participantsFoldersIndex <- 2
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
  
  # data in original space
  filesReference <- dir( pattern='sub-*' )

  # roi
  instr <- sprintf( '3dcalc -a aparc+aseg_REN_all.nii.gz -expr \u027within(a,79.95,80.05)\u027 -prefix roi_out_or.nii.gz'  ); system( instr )
  instr <- sprintf('3dresample -master %s -rmode NN -inset roi_out_or.nii.gz -prefix roi_out.nii.gz', filesReference[1] ); system( instr )

  # curvature
  instr <- 'ConvertDset -o_1D -input std.141.rh.curv.niml.dset -prepend_node_index_1D -prefix std.141.rh.curv_ppp_.1D.dset'; system( instr )
  instr <- sprintf( '3dSurf2Vol -spec std.141.AHEAD_test_rh.spec -surf_A std.141.rh.white.gii -surf_B std.141.rh.pial.gii -sv %s -grid_parent %s -map_func ave -sdata_1D std.141.rh.curv_ppp_.1D.dset -datum float -f_steps 60 -f_index voxels -prefix curvature_rh_ppp_.nii.gz', filesReference[1], filesReference[1] )
  print( instr )
  system( instr ) 
  instr <- 'ConvertDset -o_1D -input std.141.lh.curv.niml.dset -prepend_node_index_1D -prefix std.141.lh.curv_ppp_.1D.dset'; system( instr )
  instr <- sprintf( '3dSurf2Vol -spec std.141.AHEAD_test_lh.spec -surf_A std.141.lh.white.gii -surf_B std.141.lh.pial.gii -sv %s -grid_parent %s -map_func ave -sdata_1D std.141.lh.curv_ppp_.1D.dset -datum float -f_steps 60 -f_index voxels -prefix curvature_lh_ppp_.nii.gz', filesReference[1], filesReference[1] )
  print( instr )
  system( instr ) 
  # combine lh and rh curvature
  instr <- '3dcalc -a curvature_lh_ppp_.nii.gz -b curvature_rh_ppp_.nii.gz -expr \u0027a+b\u0027 -prefix curvature_both_ppp_.nii.gz'; system( instr )
  
  # thickness
  instr <- 'ConvertDset -o_1D -input std.141.rh.thickness.niml.dset -prepend_node_index_1D -prefix std.141.rh.thickness_ppp_.1D.dset'; system( instr )
  instr <- sprintf( '3dSurf2Vol -spec std.141.AHEAD_test_rh.spec -surf_A std.141.rh.white.gii -surf_B std.141.rh.pial.gii -sv %s -grid_parent %s -map_func ave -sdata_1D std.141.rh.thickness_ppp_.1D.dset -datum float -f_steps 60 -f_index voxels -prefix thickness_rh_ppp_.nii.gz', filesReference[1], filesReference[1] )
  print( instr )
  system( instr ) 
  instr <- 'ConvertDset -o_1D -input std.141.lh.thickness.niml.dset -prepend_node_index_1D -prefix std.141.lh.thickness_ppp_.1D.dset'; system( instr )
  instr <- sprintf( '3dSurf2Vol -spec std.141.AHEAD_test_lh.spec -surf_A std.141.lh.white.gii -surf_B std.141.lh.pial.gii -sv %s -grid_parent %s -map_func ave -sdata_1D std.141.lh.thickness_ppp_.1D.dset -datum float -f_steps 60 -f_index voxels -prefix thickness_lh_ppp_.nii.gz', filesReference[1], filesReference[1] )
  print( instr )
  system( instr ) 
  # combine lh and rh thickness
  instr <- '3dcalc -a thickness_lh_ppp_.nii.gz -b thickness_rh_ppp_.nii.gz -expr \u0027a+b\u0027 -prefix thickness_both_ppp_.nii.gz'; system( instr )
  
  roiFile <- read.AFNI('roi_out.nii.gz')
  depthFile <- read.AFNI('volumetricData_layering_depth.nii.gz')
  curvatureFile <- read.AFNI('curvature_both_ppp_.nii.gz')
  thicknessFile <- read.AFNI('thickness_both_ppp_.nii.gz')
  profilesFileT1 <- read.AFNI('profileData_t1w.nii.gz')
  profilesFileT2Star <- read.AFNI('profileData_t2Star.nii.gz')
  profilesFileT1Map <- read.AFNI('profileData_t1Map.nii.gz')
  profilesFileR2Star <- read.AFNI('profileData_r2star.nii.gz')
  profilesFileR1Map <- read.AFNI('profileData_r1m.nii.gz')
  profilesFileQSM <- read.AFNI('profileData_qsm.nii.gz')
    
  roiVolume <- roiFile$brk[,,,1]
  depthVolume <- depthFile$brk[,,,1]
  curvatureVolume <- curvatureFile$brk[,,,1]
  thicknessVolume <- thicknessFile$brk[,,,1]
  
  subjectName <- strsplit( filesReference[1], '_' )[[1]][1]
  
  # extract data fun
  profileFun <- function( profilesVolume, nameIn, subjectName ) {
    roiIdx <- which( depthVolume>0.4 & depthVolume<0.6 & abs(curvatureVolume) > 0.0000001 & abs(thicknessVolume) > 0.0000001 & roiVolume==1 )
    roiCoords <- coordinateFromLinearIndex( roiIdx, dim(roiVolume) )
    profilesArray <- array(999,c(length(roiIdx),dim(profilesVolume)[4]))
    for (profileIdx in 1:(dim(profilesVolume)[4])) { 
      tmpVolume <- profilesVolume[,,,profileIdx] 
      profilesArray[,profileIdx] <- tmpVolume[ roiIdx ] 
    }
    curvatureArray <- curvatureVolume[ roiIdx ]
    thicknessArray <- thicknessVolume[ roiIdx ]
    profilesArray <- cbind( curvatureArray, t(roiCoords), profilesArray )
    profilesArray <- data.frame( round( cbind( thicknessArray, profilesArray ), 5 ) )
    profilesArray <- cbind( rep( subjectName, dim(profilesArray)[1] ), rep( nameIn, dim(profilesArray)[1] ), profilesArray )
    names( profilesArray ) <- c( 'subj', 'modality', 'thickness', 'curvature','x','y','z', letters[ 1 : ( dim( profilesArray )[2] - 7 ) ] )
    return( profilesArray  )
  }
  
  profilesVolume <- profilesFileT1$brk; nameIn <- 'T1W'
  profilesT1W <- profileFun( profilesVolume, nameIn, subjectName ) 
  
  profilesVolume <- profilesFileT1Map$brk; nameIn <- 'T1Map'
  profilesT1Map <- profileFun( profilesVolume, nameIn, subjectName ) 

  profilesVolume <- profilesFileT2Star$brk; nameIn <- 'T2Star'
  profilesT2Star <- profileFun( profilesVolume, nameIn, subjectName ) 

  profilesVolume <- profilesFileQSM$brk; nameIn <- 'QSM'
  profilesQSM <- profileFun( profilesVolume, nameIn, subjectName ) 

  profilesVolume <- profilesFileR1Map$brk; nameIn <- 'R1Map'
  profilesR1Map <- profileFun( profilesVolume, nameIn, subjectName ) 

  profilesVolume <- profilesFileR2Star$brk; nameIn <- 'R2Star'
  profilesR2Star <- profileFun( profilesVolume, nameIn, subjectName ) 
  
  dfOut <- rbind( profilesT1W, profilesT1Map, profilesT2Star, profilesQSM, profilesR1Map, profilesR2Star  ) 
  
  save( dfOut, file = sprintf( '%s/%s_%s.RData', saveRoi, subjectName, appendFilename ) )
  
  # clean up
  system('rm roi_out*')
  system('rm *_ppp_*')
  
  # clean memory
  rm( list=c('profilesFileT1','profilesFileT2Star','profilesFileT1Map','profilesFileR2Star','profilesFileR1Map','profilesFileQSM') ); 
  rm( list=c('profilesVolume') ); 
  gc()
  
}