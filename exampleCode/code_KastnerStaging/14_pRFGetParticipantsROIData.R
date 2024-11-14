rm(list=ls())

# toolbox and functions
generalPurposeDir <- Sys.getenv( x='AFNI_TOOLBOXDIRGENERALPURPOSE' )
afniInstallDir <- Sys.getenv( x='AFNI_INSTALLDIR' ) 
source( sprintf('%s/scaleData.R', generalPurposeDir) )
source( sprintf('%s/AFNIio.R', afniInstallDir ) )
library( pracma )
library( abind )
library( circular )

# folders 
mainFolder <- '/analyse/Project0226/KastnerModel/staging_area_Kastner'
anatomyDir <- sprintf( '%s/anatomies_KastnerClassic_Freesurfer', mainFolder )
modelOutputFolder <- sprintf('%s/modelsOutput/pRF', mainFolder )
surfaceFolder <- sprintf('%s/anatomies_KastnerClassic_Freesurfer/surfaceAtlases/suma_MNI152_2009', mainFolder)
roiSurfaceFolderMNI <- sprintf('%s/probatlas_v4/ProbAtlas_v4/subj_surf_all', surfaceFolder)
resultsFolder <- sprintf('%s/results', mainFolder ) 

setwd( modelOutputFolder )
print( sprintf('current folder: %s', getwd() ) )
modelParticipants <- dir( modelOutputFolder )
print( sprintf('model participants:' ) )
modelParticipants

setwd( anatomyDir )
print( sprintf('current folder: %s', getwd() ) )
anatomyFolders <- dir( anatomyDir, pattern='*_ANATOMY' )
print( sprintf('anatomy folders:' ) )
anatomyFolders

# look for corresponding anatomy folders
anatomyFolders_id <- rep('a',length(anatomyFolders))
for (nAnat in 1:length(anatomyFolders)) {
  anatomyFolders_id[nAnat] <- strsplit( anatomyFolders, '_' )[[nAnat]][1]
}
anatomyFolders_matching <- rep(0,length(modelParticipants))
for ( nPart in 1:length( modelParticipants ) ) {
  partId <- strsplit( modelParticipants[ nPart ], '_' )[[1]][1]
  anatIdx <- which( anatomyFolders_id == partId )
  print( anatIdx )
  anatomyFolders_matching[nPart] <- anatIdx
}

# check the alignment
selectedAnatFolders <- anatomyFolders[ anatomyFolders_matching ]
data.frame( modelParticipants, selectedAnatFolders )

runCodeFlag <- 1

Sys.setenv(OMP_NUM_THREADS='8')
Sys.getenv('OMP_NUM_THREADS')

#function to load models
loadModel <- function( modelFileName, roiLeft, roiRight, participant ) {
  
  roiNames <- c('V1v','V1d','V2v','V2d','V3v','V3d','hV4','VO1','VO2','PHC1','PHC2','MST','hMT',	   
                'LO2','LO1','V3b','V3a','IPS0','IPS1','IPS2','IPS3','IPS4','IPS5','SPL1','FEF')	
  
  # participant <- 'ASM16'
  # modelFileName <- allParams[1]
  modelFile <- read.AFNI( filename = modelFileName )
  modelVolume <- modelFile$brk
  dfOut <- c()
  
  for ( nHemi in c('lh','rh') ) { # nHemi <- 'lh'
    
    if ( nHemi=='lh') { roiFile <- roiLeft; roiIdx <- 1 }
    if ( nHemi=='rh') { roiFile <- roiRight; roiIdx <- 2 }
    roiVolume <- roiFile$brk
    
    testDimensionFlag <- sum( dim( roiVolume )[c(1,2,3)] == dim( modelVolume )[c(1,2,3)] )
    if ( testDimensionFlag != 3 ) { break; print('dimensions mismatch') }
    
    for ( nRoi in 1:25 ) { #25 roi indexes nRoi <- 1
      idxRoi <- which( roiVolume == nRoi )
      dfRoi <- c()
      for (nBrik in 1:dim(modelVolume)[4]) { #nBrik <- 2
        print( sprintf('loading model: %s, brik: %d, ROI: %d, hemisphere: %s; hemiIdx: %d', modelFileName, nBrik, nRoi, nHemi, roiIdx ) )
        modelVolume_individualBrik <- modelVolume[,,,nBrik]
        dfRoi <- cbind( dfRoi, modelVolume_individualBrik[ idxRoi ] )
      }
      dfRoi <- cbind( rep( roiIdx, dim( dfRoi )[1] ), rep( nRoi, dim( dfRoi )[1] ), dfRoi )
      dfOut <- rbind( dfOut, dfRoi )
    }
  }
  dfOut <- round( dfOut, 4 )
  dfOutDF <- data.frame( rep( participant, dim(dfOut)[1] ),
                         rep( c('lh','rh'), array( table( dfOut[,1] ) ) ),
                         dfOut )
  modelString <- gsub( pattern='\"', '', modelFile$NI_head$BRICK_LABS$dat )
  modelString01 <- strsplit( modelString, ' ' )[[1]][2]
  modelString02 <- strsplit( modelString01, '~' )[[1]]
  dfNames <- c( 'participant', 'hemiName', 'hemiIdx', 'roiIdx', modelString02 )  
  names( dfOutDF ) <- dfNames
  dfOutDF$ROI_name <- roiNames[ dfOutDF$roiIdx ]
  
  nameTmp <- dfOutDF$ROI_name
  dfOutDF$ROI_name_basic <- ifelse( or( nameTmp=='V1v', nameTmp=='V1d' ), 'V1',
                                    ifelse( or( nameTmp=='V2v', nameTmp=='V2d' ), 'V2', 
                                            ifelse( or( nameTmp=='V3v', nameTmp=='V3d' ), 'V3', nameTmp 
                                            ) ) ) 
  print( 'roi name contingency table, unsorted:' )
  print( table( dfOutDF$ROI_name_basic, dfOutDF$ROI_name ) )
  
  dfOutDF$ROI_name <- factor( dfOutDF$ROI_name, levels=unique( dfOutDF$ROI_name ) ) #V1v, V1d, V2v, V2d etc etc (Wang2015 atlas names)
  dfOutDF$ROI_name_basic <- factor( dfOutDF$ROI_name_basic, levels=unique( dfOutDF$ROI_name_basic ) ) #V1, V2, V3 etc etc (Wang2015 atlas but without distinguishing between ventral and dorsal early visual cortex)
  
  print( 'roi name contingency table, sorted:' )
  print( table( dfOutDF$ROI_name_basic, dfOutDF$ROI_name ) ) # contingency table to check roi naming
  
  return( dfOutDF )
}

# move to the folder of the first participant, where ALL modelling results are stored (also individual Cw and CCw models are stored)
# just to get model names
setwd( modelOutputFolder )
setwd( modelParticipants[ 1 ] )
print('########')
print( getwd() )
print('########')

paramsFilepRF <- list.files( path=getwd(), pattern='*params*', recursive = FALSE, full.names = FALSE  )

allParamsTemp <- c( paramsFilepRF )
allParams <- c()
for ( n in 1:length(allParamsTemp) ) {
  tempStr <- strsplit(allParamsTemp[n],'AFH28')[[1]][2]
  allParams <- c( allParams, tempStr )
}

#### create output folder ####
setwd( resultsFolder )
dirToCheck <- sprintf('pRF' )
if ( dir.exists( dirToCheck ) ) {
  instr <- sprintf( 'rm -R %s', dirToCheck ) 
  print( instr )
  if (runCodeFlag==1) { system( instr ) }  
}  
if ( !dir.exists( dirToCheck ) ) {
  instr <- sprintf('mkdir %s', dirToCheck)
  print( instr )
  if (runCodeFlag==1) { system( instr ) }  
}         


for ( nModels in 1:length( allParams ) ) { # nModels <- 1
  outModelDf <- c()
  for ( nPart in 1:length( modelParticipants )  ) { # nPart <- 1 # length( modelParticipants )
    
    # move to the right anatomy folder, where modelling results are stored
    setwd( anatomyDir )
    setwd( selectedAnatFolders[ nPart ] )
    setwd( 'FreeSeg_result/SUMA' )
    setwd( 'roisWang2015_pRF' ) # use rois from pRF data, however all the EPIS are coregistered to the same anatomy, so we can compare between models
    print('########')
    print( getwd() )
    print('########')
    
    # get roi data:
    roiLeft <- read.AFNI( filename='WangRoi2015.lh.nii.gz' )
    roiRight <- read.AFNI( filename='WangRoi2015.rh.nii.gz' )
    
    # move to the folder individual participant, where ALL modelling results are stored (also individual Cw and CCw models are stored)
    setwd( modelOutputFolder )
    setwd( modelParticipants[ nPart ] )
    print('########')
    print( getwd() )
    print('########')
    
    participant <- modelParticipants[ nPart ]
    modelFileName <- paste( participant, allParams[ nModels ], sep='')
    outModelDfTemp <- loadModel( modelFileName, roiLeft, roiRight, participant )
    outModelDfTemp$modelName <- allParams[ nModels ]
    outModelDf <- rbind( outModelDf, outModelDfTemp )
    
  }
  save( outModelDf, file = sprintf( '%s/%s/%s.RData', resultsFolder, 'pRF', allParams[ nModels ] ) )
}

