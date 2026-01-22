rm( list=ls() )
source('~/abin/AFNIio.R')

debugFlag  <- 0
if (debugFlag==1) {
  mainFolder <- '/media/alessiofracasso/DATADRIVE1/Flanker'
  inputFolder <- 'derivatives/processing_afni_denoised'
  outputFolder <- 'derivatives/resultsGlasser'
  outputFileName <- 'results_Glasser_ORIG_Denoised.RData'
  roiFilenameLH_input <- 'lh.Glasser_HCP_epiResampled_' 
  roiFilenameRH_input <- 'rh.Glasser_HCP_epiResampled_' 
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
      roiVolume <- roiFile$brk
      roiNames <-  c('V1','MST','V6','V2','V3','V4','V8','4','3b','FEF','PEF','55b','V3A','RSC','POS2','V7','IPS1','FFC','V3B','LO1','LO2','PIT','MT','A1','PSL','SFL',
                     'PCV','STV','7Pm','7m','POS1','23d','v23ab','d23ab','31pv','5m','5mv','23c','5L','24dd','24dv','7AL','SCEF','6ma','7Am','7Pl','7PC','LIPv','VIP','MIP',
                     '1','2','3a','6d','6mp','6v','p24pr','33pr','a24pr','p32pr','a24','d32','8BM','p32','10r','47m','8Av','8Ad','9m','8BL','9p','10d','8C','44','45','47l',
                     'a47r','6r','IFJa','IFJp','IFSp','IFSa','p9-46v','46','a9-46v','9-46d','9a','10v','a10p','10pp','11l','13l','OFC','47s','LIPd','6a','i6-8','s6-8','43',
                     'OP4','OP1','OP2-3','52','RI','PFcm','PoI2','TA2','FOP4','MI','Pir','AVI','AAIC','FOP1','FOP3','FOP2','PFt','AIP','EC','PreS','H','ProS','PeEc','STGa',
                     'PBelt','A5','PHA1','PHA3','STSda','STSdp','STSvp','TGd','TE1a','TE1p','TE2a','TF','TE2p','PHT','PH','TPOJ1','TPOJ2','TPOJ3','DVT','PGp','IP2','IP1','IP0',
                     'PFop','PF','PFm','PGi','PGs','V6A','VMV1','VMV3','PHA2','V4t','FST','V3CD','LO3','VMV2','31pd','31a','VVC','25','s32','pOFC','PoI1','Ig','FOP5','p10p','p47r',
                     'TGv','MBelt','LBelt','A4','STSva','TE1m','PI','a32pr','p24')
      
      roiLabelsIdx00 <- sort( unique( array( roiVolume ) ) ) # sort because roi indexes MUST be in register with roi names
      roiLabelsIdx01 <- roiLabelsIdx00[ roiLabelsIdx00 != 0 ] # -1 to remove 0 (no ROI)
      nRois <- length( roiLabelsIdx01  )
      
      #dataPartSideTemp <- data.frame( array( 999, c( length(statNames)*length(roiNames), 1 ) ) )
      #names( dataPartSideTemp ) <- c('medianValue')
      #dataPartSideTemp$side <- side
      #dataPartSideTemp$participantID <- partNameClean
      #dataPartSideTemp$roi <- 'tempRoi'
      #dataPartSideTemp$coefficient <- 'tempCoefficient'
      #dataPartSideTemp$roiIdx <- 999
      
      dataPartSideTemp <- data.frame( voxelValue=double(),
      				      voxelIdx=double(),
                                      side=factor(),
                                      participantID=factor(),
                                      roi=factor(),
                                      coefficient=double(),
                                      roiIdx=double() )
      
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
          if (length(roiValuesTemp_check)==0) { roiValuesTemp_check <- 0 }
          
          voxelIdx_roi <- seq( 1, length(roiValuesTemp_check) )
          
          computedMedian <- round( median( roiValuesTemp_check ), 4 )
          
          print( computedMedian )
          print( sprintf( 'length(roiValuesTemp_check): %d', length(roiValuesTemp_check) ) )
          print( sprintf( 'length( voxelIdx_roi ): %d', length( voxelIdx_roi ) ) )
          
          #dataPartSideTemp$medianValue[ lineCounter ] <- computedMedian
          #dataPartSideTemp$roi[ lineCounter ] <- roiNames[ roiCounter ]
          #dataPartSideTemp$coefficient[ lineCounter ] <- statNames[ coeffCounter ]
          #dataPartSideTemp$roiIdx[ lineCounter ] <- roiCounter
          
          dataPartSideTempLoop <- data.frame( voxelValue=roiValuesTemp_check,
          					voxelIdx=voxelIdx_roi,
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
inputFolder <- 'derivatives/processing_afni_denoised'
outputFolder <- 'derivatives/resultsGlasser'
outputFileName <- 'results_Glasser_ORIG_Denoised_singleVoxel.RData'
roiFilenameLH_input <- 'lh.Glasser_HCP_epiResampled_' 
roiFilenameRH_input <- 'rh.Glasser_HCP_epiResampled_' 
MNI_flag <- 0
getResultsGlasser( mainFolder, inputFolder, outputFolder, outputFileName, roiFilenameLH_input, roiFilenameRH_input, MNI_flag )


#fileToLoad <- sprintf('%s/%s/results_Glasser.RData', mainDir, outputDir)
#load( file=fileToLoad )

#table( outputDataset$participantID, outputDataset$roi )
#table( outputDataset$participantID, outputDataset$roi, outputDataset$side )
