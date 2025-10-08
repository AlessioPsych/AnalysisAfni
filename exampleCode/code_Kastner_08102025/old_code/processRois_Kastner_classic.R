rm( list=ls() );
mainDir <- '/analyse/Project0226/dataRepository'
targetDir_pRF <- '/analyse/Project0226/KastnerModel/data_KastnerClassic_pRF'
targetDir_KastnerClassic <- '/analyse/Project0226/KastnerModel/data_KastnerClassic'
anatomiesDir_whole <- '/analyse/Project0226/KastnerModel/anatomies_KastnerClassic'
setwd( mainDir )

pRFDir <- dir( targetDir_pRF )
KCDir <- dir( targetDir_KastnerClassic )
anatomiesDir <- dir( anatomiesDir_whole )

data.frame( pRFDir, KCDir, anatomiesDir )

generalPurposeDir <- Sys.getenv( x='AFNI_TOOLBOXDIRGENERALPURPOSE' )
afniInstallDir <- Sys.getenv( x='AFNI_INSTALLDIR' ) 
source( sprintf('%s/scaleData.R', generalPurposeDir) )
source( sprintf('%s/AFNIio.R', afniInstallDir ) )
library( pracma )
library( abind )
library( circular )

for (nSubjects in 1:length(anatomiesDir)) { #subjId <- 1
  
  subjId <- nSubjects
  
  # clean data on the Kastner classic pRF folder
  setwd( targetDir_pRF ); getwd()
  setwd( pRFDir[ subjId ] ); getwd()
  system( 'rm -R roi_epi/' )
  system( 'mkdir roi_epi' )
  
  # clean data on the Kastner classic folder
  setwd( targetDir_KastnerClassic ); getwd()
  setwd( KCDir[ subjId ] ); getwd()
  system( 'rm -R roi_epi/' )
  system( 'mkdir roi_epi' )
  
  # clean data on the anatomy folder
  setwd( anatomiesDir_whole )
  setwd( anatomiesDir[ subjId ] ); getwd()
  system( 'rm -R roi_volume/' )
  system( 'mkdir roi_volume' )

  # convert rois from surface to volume
  setwd('roi')
  roiFilesAll <- dir( pattern='*.1D.roi' )
  keepRois <- function (idx, roiFilesAll) { is.null( strfind( 'borders', strsplit( roiFilesAll[idx], '_' )[[1]][1] ) ) }  
  keepRoisFlag <- sapply( 1:length( roiFilesAll ), keepRois, roiFilesAll=roiFilesAll  )
  roiFiles <- roiFilesAll[ keepRoisFlag ]
  print('..........................')
  print('..........................')
  print('rois to process...........')
  print('..........................')
  print('..........................')
  print( roiFiles )
  setwd( '..' )
  
  for ( nRois in 1:length(roiFiles)) { #length(roiFiles) #nRois <- 1

    print('.............................')
    print('.............................')
    print('.............................')
    print( sprintf('processing roi: %s, participant: %s, participant (pRF): %s; ', roiFiles[nRois], KCDir[ subjId ], pRFDir[ subjId ] ) )
    print('.............................')
    print('.............................')
    print('.............................')
    
    instr <- sprintf('cp roi/%s %s', roiFiles[nRois], roiFiles[nRois]  ); system( instr )
    side <- strsplit( roiFiles[nRois], '_' )[[1]][2]
    if (side=='left') {
      instr <- sprintf('surf2vol.sh %s boundary01 volumetricData_layering_depth.nii.gz 30 surfaces_folder_left/', roiFiles[nRois] ); system( instr )
    }
    if (side=='right') {
      instr <- sprintf('surf2vol.sh %s boundary01 volumetricData_layering_depth.nii.gz 30 surfaces_folder_right/', roiFiles[nRois] ); system( instr )
    }
    roiFilesOutcome <- dir( pattern=strsplit( roiFiles[nRois], '[.]' )[[1]][1]  )
    instr <- sprintf('gunzip %s', roiFilesOutcome[2] ); system( instr )
    roiFilesOutcome <- dir( pattern=strsplit( roiFiles[nRois], '[.]' )[[1]][1]  )
    for (nFiles in 1:length(roiFilesOutcome) ) {
      instr <- sprintf( 'mv %s roi_volume/%s', roiFilesOutcome[nFiles], roiFilesOutcome[nFiles] ); system( instr )
    }
    setwd('roi_volume')
    instr <- sprintf( '3dAllineate -source %s -prefix %s/%s/roi_epi/%s -1Dmatrix_apply %s/%s/coregistration_kastner/invMat.1D -master %s/%s/meanEpi4Coreg_kastner.nii.gz -final NN', roiFilesOutcome[2], targetDir_KastnerClassic, KCDir[ subjId ], roiFilesOutcome[2], targetDir_KastnerClassic, KCDir[ subjId ], targetDir_KastnerClassic, KCDir[ subjId ]  )
    system( instr )
    instr <- sprintf( '3dAllineate -source %s -prefix %s/%s/roi_epi/%s -1Dmatrix_apply %s/%s/coregistration_bars/invMat.1D -master %s/%s/meanEpi4Coreg_bars.nii.gz -final NN', roiFilesOutcome[2], targetDir_pRF, pRFDir[ subjId ], roiFilesOutcome[2], targetDir_pRF, pRFDir[ subjId ], targetDir_pRF, pRFDir[ subjId ]  )
    system( instr )
    setwd('..')
  }
  
}


