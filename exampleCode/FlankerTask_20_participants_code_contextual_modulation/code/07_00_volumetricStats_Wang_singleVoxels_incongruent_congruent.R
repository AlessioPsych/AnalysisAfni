rm( list=ls() )
source('~/abin/AFNIio.R')

debugFlag  <- 0
if (debugFlag==1) {
  mainFolder <- '/home/fracasso/Data/openNeuro/ds000102'
  inputFolder <- 'derivatives/processing_afni_denoised_incongruent_congruent'
  outputFolder <- 'derivatives/resultsWang'
  outputFileName <- 'results_Wang_ORIG_Denoised.RData'
  roiFilenameLH_input <- 'lh.Wang_2015_epiResampled_' 
  roiFilenameRH_input <- 'rh.Wang_2015_epiResampled_' 
  MNI_flag <- 0
}

getResultsGlasser <- function( mainFolder, inputFolder, outputFolder, outputFileName, roiFilenameLH_input, roiFilenameRH_input, MNI_flag ) {
  
  mainDir <- mainFolder
  dataDir <- inputFolder
  outputDir <- outputFolder
  setwd( mainDir )
  setwd( dataDir )
  
  partDirs <- list.dirs(full.names = FALSE, recursive = FALSE)
  partDirs <- partDirs[grepl("^*sub-", partDirs)]
  print( 'participant dirs...' )
  print(partDirs)
  
  statFilename <- 'stats.'
  tsnrFilename <- 'TSNR.'
  roiFilenameLH <- roiFilenameLH_input
  roiFilenameRH <- roiFilenameRH_input
  #  roiFilenameLH <- 'lh.Glasser_HCP_epiResampled_' 
  #  roiFilenameRH <- 'rh.Glasser_HCP_epiResampled_' 
  
  #flagOutput <- 1
  outputDataset <- data.frame( voxelValue=double(), 
  			       voxelIdx=double(),
                               side=factor(),
                               participantID=factor(),
                               roi=factor(),
                               coefficient=factor(),
                               roiIdx=double())
  for ( i in  1:length( partDirs ) ) { # i <- 1; side <- 'lh' length( partDirs )
    for ( side in c('lh','rh') ) {
      setwd( mainDir )
      setwd( dataDir )
      setwd( partDirs[ i ] )
      print( getwd() )
      print( sprintf('side %s', side) )

      partNameClean <- strsplit( partDirs[ i ], '.results' )[[1]][1]

      if ( side=='lh') { roiFilenameTemp <- roiFilenameLH }
      if ( side=='rh') { roiFilenameTemp <- roiFilenameRH }
      roiFilenameLoop <- sprintf( '%s%s.nii.gz', roiFilenameTemp, partNameClean )
      
      if (MNI_flag==0) { statFilenameLoop <- sprintf( '%s%s+orig', statFilename, partNameClean ) }
      if (MNI_flag==1) { statFilenameLoop <- sprintf( '%s%s+tlrc', statFilename, partNameClean ) }
      if (MNI_flag==0) { tsnrFilenameLoop <- sprintf( '%s%s+orig', tsnrFilename, partNameClean ) }
      if (MNI_flag==1) { tsnrFilenameLoop <- sprintf( '%s%s+tlrc', tsnrFilename, partNameClean ) }
      if (MNI_flag==0) { retinoFilenameLoop <- sprintf( '%s.Benson_Retinotopy_epiResampled_%s.nii.gz', side, partNameClean ) }
      
      # delete concatenate stats and tsnr in a single volume, if it exists
      if ( file.exists( 'concatenatedDelMe.nii.gz' ) ) {
        instr <- sprintf( 'rm concatenatedDelMe.nii.gz' )
        system( instr )
      }
      
      # concatenate stats, tsnr and retino in a single volume
      instr <- sprintf( '3dTcat -prefix concatenatedDelMe.nii.gz %s %s %s', statFilenameLoop, tsnrFilenameLoop, retinoFilenameLoop )
      system( instr )
      
      #if ( side=='lh') { roiFilenameTemp <- roiFilenameLH }
      #if ( side=='rh') { roiFilenameTemp <- roiFilenameRH }
      #roiFilenameLoop <- sprintf( '%s%s.nii.gz', roiFilenameTemp, partNameClean )
      
      print( statFilenameLoop ) 
      print( roiFilenameLoop ) 
      print( retinoFilenameLoop ) 
      
      statFile_names <- read.AFNI( statFilenameLoop )
      statNames <-  strsplit( statFile_names$NI_head$BRICK_LABS$dat, '~' )[[1]]
      statNames[ length( statNames ) + 1 ] <- 'tsnr'
      statNames[ length( statNames ) + 1 ] <- 'voxelIndx'
      statNames[ length( statNames ) + 1 ] <- 'R2'
      statNames[ length( statNames ) + 1 ] <- 'rfSize'
      statNames[ length( statNames ) + 1 ] <- 'meanVol'
      statNames[ length( statNames ) + 1 ] <- 'gain'
      statNames[ length( statNames ) + 1 ] <- 'ecc'
      statNames[ length( statNames ) + 1 ] <- 'ang'      
      
      statFile <- read.AFNI( 'concatenatedDelMe.nii.gz' ) #with tsnr concatenated
      statVolume <- statFile$brk
      
      # delete concatenate stats and tsnr in a single volume, if it exists
      if ( file.exists( 'concatenatedDelMe.nii.gz' ) ) {
        instr <- sprintf( 'rm concatenatedDelMe.nii.gz' )
        print( instr )
        system( instr )
      }
      
      roiFile <- read.AFNI( roiFilenameLoop )
      roiVolume <- roiFile$brk[,,,2] # wang atlas is in the second volume
      roiNames <-  c('V1v',	'V1d', 'V2v', 'V2d', 'V3v', 'V3d', 'hV4', 'VO1', 'VO2', 'PHC1', 'PHC2',	    
      'MST', 'hMT', 'LO2', 'LO1', 'V3b', 'V3a', 'IPS0', 'IPS1', 'IPS2', 'IPS3', 'IPS4', 'IPS5', 'SPL1','FEF')	
      
      roiLabelsIdx00 <- sort( unique( array( roiVolume ) ) ) # sort because roi indexes MUST be in register with roi names
      roiLabelsIdx01 <- roiLabelsIdx00[ roiLabelsIdx00 != 0 ] # -1 to remove 0 (no ROI)
      nRois <- length( roiLabelsIdx01  )
      
      dataPartSideTemp <- data.frame( voxelValue=double(),
      				   voxelIdx=double(), 
                                   side=factor(),
                                   participantID=factor(),
                                   roi=factor(),
                                   coefficient=factor(),
                                   roiIdx=double())
      
      #lineCounter <- 1
      for ( roiCounter in 1:nRois) { # roiCounter <- 1
        for ( coeffCounter in 1:dim(statVolume)[4] ) { # coeffCounter <- 1
          print( partNameClean )
          print( sprintf( 'roi processed: %d, coeff processed: %d', roiCounter, coeffCounter ) )
          print( sprintf('side %s', side) )
          idxRoiTemp <- which( roiVolume == roiCounter )
          statVolumeTemp <- statVolume[,,,coeffCounter]
          roiValuesTemp_check <- statVolumeTemp[ idxRoiTemp ]
          
          #roiValuesTemp_check <- roiValuesTemp[ abs(roiValuesTemp) > 0.00001 ] # keep only values effectively larger (or smaller) than zero
          if ( length(roiValuesTemp_check)==0 & length(idxRoiTemp)==0 ) { 
          	roiValuesTemp_check <- 0 
          	idxRoiTemp <- 0
          }
          
          voxelIdx_roi <- seq( 1, length(roiValuesTemp_check) )
                    
          computedMedian <- round( median( roiValuesTemp_check ), 4 )
          
          print( computedMedian )
          print( sprintf( 'length(roiValuesTemp_check): %d', length(roiValuesTemp_check) ) )
          print( sprintf( 'length( voxelIdx_roi ): %d', length( voxelIdx_roi ) ) )
          
          dataPartSideTempLoop <- data.frame( voxelValue=roiValuesTemp_check,
          				      voxelIdx=idxRoiTemp,	
                                              side=side,
                                              participantID=partNameClean,
                                              roi=roiNames[ roiCounter ],
                                              coefficient=statNames[ coeffCounter ],
                                              roiIdx=roiCounter )
          
          dataPartSideTemp <- rbind( dataPartSideTemp, dataPartSideTempLoop )
          
          #lineCounter <- lineCounter + 1
        }
      }
      
      outputDataset <- rbind( outputDataset, dataPartSideTemp )
      #if ( flagOutput==0 ) { outputDataset <- rbind( outputDataset, dataPartSideTemp ) }
      #if ( flagOutput==1 ) { outputDataset <- dataPartSideTemp; flagOutput <- 0 }
      
    }
  }
  
  fileToSave <- sprintf('%s/%s/%s', mainDir, outputDir, outputFileName )
  save( outputDataset, file=fileToSave )
  
}

mainFolder <- '/home/fracasso/Data/openNeuro/ds000102'
inputFolder <- 'derivatives/processing_afni_denoised_incongruent_congruent'
outputFolder <- 'derivatives/resultsWang'
outputFileName <- 'results_Wang_ORIG_Denoised_singleVoxel_incongruent_congruent.RData'
roiFilenameLH_input <- 'lh.Wang_2015_epiResampled_' 
roiFilenameRH_input <- 'rh.Wang_2015_epiResampled_' 
MNI_flag <- 0
getResultsGlasser( mainFolder, inputFolder, outputFolder, outputFileName, roiFilenameLH_input, roiFilenameRH_input, MNI_flag )

#mainFolder <- '/media/alessiofracasso/DATADRIVE1/Flanker'
#inputFolder <- 'derivatives/processing_afni_MNI_denoised'
#outputFolder <- 'derivatives/resultsWang'
#outputFileName <- 'results_Wang_MNI_Denoised.RData'
#roiFilenameLH_input <- 'lh.Wang_2015_MNI_epiResampled_' 
#roiFilenameRH_input <- 'rh.Wang_2015_MNI_epiResampled_' 
#MNI_flag <- 1
#getResultsGlasser( mainFolder, inputFolder, outputFolder, outputFileName, roiFilenameLH_input, roiFilenameRH_input, MNI_flag )

#mainFolder <- '/media/alessiofracasso/DATADRIVE1/Flanker'
#inputFolder <- 'derivatives/processing_afni_no_denoised'
#outputFolder <- 'derivatives/resultsWang'
#outputFileName <- 'results_Wang_ORIG_No_Denoised.RData'
#roiFilenameLH_input <- 'lh.Wang_2015_epiResampled_' 
#roiFilenameRH_input <- 'rh.Wang_2015_epiResampled_' 
#MNI_flag <- 0
#getResultsGlasser( mainFolder, inputFolder, outputFolder, outputFileName, roiFilenameLH_input, roiFilenameRH_input, MNI_flag )

#mainFolder <- '/home/fracasso/Data/openNeuro/ds000102'
#inputFolder <- 'derivatives/processing_afni_MNI_no_denoised'
#outputFolder <- 'derivatives/resultsWang'
#outputFileName <- 'results_Wang_MNI_No_Denoised.RData'
#roiFilenameLH_input <- 'lh.Wang_2015_MNI_epiResampled_' 
#roiFilenameRH_input <- 'rh.Wang_2015_MNI_epiResampled_' 
#MNI_flag <- 1
#getResultsGlasser( mainFolder, inputFolder, outputFolder, outputFileName, roiFilenameLH_input, roiFilenameRH_input, MNI_flag )

#fileToLoad <- sprintf('%s/%s/results_Glasser.RData', mainDir, outputDir)
#load( file=fileToLoad )

#table( outputDataset$participantID, outputDataset$roi )
#table( outputDataset$participantID, outputDataset$roi, outputDataset$side )
