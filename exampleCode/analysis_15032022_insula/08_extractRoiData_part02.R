rm( list=ls() )
afniDir <- Sys.getenv(x='AFNI_INSTALLDIR')
afniSurfaces <- Sys.getenv(x='AFNI_TOOLBOXDIRSURFACES')
source( sprintf( '%s/AFNIio.R', afniDir ) )
source( sprintf( '%s/coordinateFromLinearIndex.r', afniSurfaces ) )

saveRoiFolder <- '/analyse/Project0165/Fracasso_Anatomy/AHEAD_database/10007840/results/insulaSurface'

mainDir <- '/analyse/Project0165/Fracasso_Anatomy/AHEAD_database/10007840'
partDir <- 'part2/'

extractSurfaceRoiData <- function( modality, roiIdx, hemi ) {
  
  if (hemi=='LH') {
    atlas <- read.table('std.141.aparc+aseg_REN_all_midpoint.lh.1D.dset', header=FALSE)
    names( atlas ) <- c('node', 'index', 'i', 'j', 'k', 'vals', 'atlasIdx' )
    
    curvature <- read.table('std.141.lh.curv_ppp_.1D.dset', header=FALSE)
    names( curvature ) <- c('node', 'curvature' )
    
    sulc <- read.table('std.141.lh.sulc_ppp_.1D.dset', header=FALSE)
    names( sulc ) <- c('node', 'sulc' )
    
    thickness <- read.table('std.141.lh.thickness_ppp_.1D.dset', header=FALSE)
    names( thickness ) <- c('node', 'thickness' )
    
    profiles_surf <- read.table( sprintf( 'std.141.%s_segvals_10.lh.1D.dset', modality ), header=FALSE )
    
    names( profiles_surf ) <- c('node', 'index', 'i', 'j', 'k', 'vals', 'V0', 'V1', 'V2', 'V3', 'V4', 'V5', 'V6', 'V7', 'V8', 'V9' )
    str( profiles_surf )
    
    dfOutSel <- data.frame( profiles_surf[ atlas$atlasIdx == roiIdx, ],
                            curvature$curvature[ atlas$atlasIdx == roiIdx ], 
                            thickness$thickness[ atlas$atlasIdx == roiIdx ],
                            sulc$sulc[ atlas$atlasIdx == roiIdx ],
                            rep(hemi, sum( atlas$atlasIdx == roiIdx ) ),
                            rep(modality, sum( atlas$atlasIdx == roiIdx ) ) )
    
    names( dfOutSel ) <- c('node', 'index', 'i', 'j', 'k', 'vals', 'V0', 'V1', 'V2', 'V3', 'V4', 'V5', 'V6', 'V7', 'V8', 'V9', 'curvature', 'thickness', 'sulc', 'hemi', 'modality' )
    str( dfOutSel )
  }
  if (hemi=='RH') {
    atlas <- read.table('std.141.aparc+aseg_REN_all_midpoint.rh.1D.dset', header=FALSE)
    names( atlas ) <- c('node', 'index', 'i', 'j', 'k', 'vals', 'atlasIdx' )
    
    curvature <- read.table('std.141.rh.curv_ppp_.1D.dset', header=FALSE)
    names( curvature ) <- c('node', 'curvature' )
    
    sulc <- read.table('std.141.rh.sulc_ppp_.1D.dset', header=FALSE)
    names( sulc ) <- c('node', 'sulc' )
    
    thickness <- read.table('std.141.rh.thickness_ppp_.1D.dset', header=FALSE)
    names( thickness ) <- c('node', 'thickness' )
    
    profiles_surf <- read.table( sprintf( 'std.141.%s_segvals_10.rh.1D.dset', modality ), header=FALSE )
    
    names( profiles_surf ) <- c('node', 'index', 'i', 'j', 'k', 'vals', 'V0', 'V1', 'V2', 'V3', 'V4', 'V5', 'V6', 'V7', 'V8', 'V9' )
    str( profiles_surf )
    
    dfOutSel <- data.frame( profiles_surf[ atlas$atlasIdx == roiIdx, ],
                            curvature$curvature[ atlas$atlasIdx == roiIdx ], 
                            thickness$thickness[ atlas$atlasIdx == roiIdx ],
                            sulc$sulc[ atlas$atlasIdx == roiIdx ], 
                            rep(hemi, sum( atlas$atlasIdx == roiIdx ) ),
                            rep(modality, sum( atlas$atlasIdx == roiIdx ) ) )
    
    names( dfOutSel ) <- c('node', 'index', 'i', 'j', 'k', 'vals', 'V0', 'V1', 'V2', 'V3', 'V4', 'V5', 'V6', 'V7', 'V8', 'V9', 'curvature', 'thickness', 'sulc', 'hemi', 'modality' )
    str( dfOutSel )
  }
  
  return( dfOutSel )
  
}

setwd( mainDir )
setwd( partDir )
getwd()
participantsFolders <- dir()
for ( participantsFoldersIndex in seq( 1,length( participantsFolders ) ) ) { #seq( 1,length( participantsFolders )
  
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
  getwd()
  
  lhT1w <- extractSurfaceRoiData( 't1w', 80, 'LH' )
  rhT1w <- extractSurfaceRoiData( 't1w', 115, 'RH' )
  lhT1map <- extractSurfaceRoiData( 't1map', 80, 'LH' )
  rhT1map <- extractSurfaceRoiData( 't1map', 115, 'RH' )
  lhR1map <- extractSurfaceRoiData( 'r1map', 80, 'LH' )
  rhR1map <- extractSurfaceRoiData( 'r1map', 115, 'RH' )
  lhQSM <- extractSurfaceRoiData( 'qsm', 80, 'LH' )
  rhQSM <- extractSurfaceRoiData( 'qsm', 115, 'RH' )
  lhT2star <- extractSurfaceRoiData( 't2starmap', 80, 'LH' )
  rhT2star <- extractSurfaceRoiData( 't2starmap', 115, 'RH' )
  lhR2star <- extractSurfaceRoiData( 'r2starmap', 80, 'LH' )
  rhR2star <- extractSurfaceRoiData( 'r2starmap', 115, 'RH' )
  
  dataOut <- rbind( lhT1w, rhT1w, lhT1map, rhT1map, lhR1map, rhR1map, lhQSM, rhQSM, lhT2star, rhT2star, lhR2star, rhR2star )
  
  str(dataOut)
  
  save( dataOut, file = sprintf( '%s/%s.RData', saveRoiFolder, participantsFolders[ participantsFoldersIndex ] ) )
  
}