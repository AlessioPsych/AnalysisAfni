rm(list=ls())
mainFolder <- '/analyse/Project0226/KastnerModel/staging_area_Kastner'
epiDir <- sprintf( '%s/epi_data', mainFolder )
epiDir_extraData <- sprintf( '%s/epi_data_extra', mainFolder )
anatomyDir <- sprintf( '%s/anatomies_KastnerClassic_Freesurfer', mainFolder )

setwd( epiDir )
print( sprintf('current folder: %s', getwd() ) )
epiSessions_orig <- dir( epiDir )
print( sprintf('epi sessions:' ) )
epiSessions_orig

setwd( epiDir_extraData )
print( sprintf('current folder: %s', getwd() ) )
epiSessions_extra <- dir( epiDir_extraData )
print( sprintf('epi sessions extra:' ) )
epiSessions_extra

setwd( anatomyDir )
print( sprintf('current folder: %s', getwd() ) )
anatomyFolders <- dir( anatomyDir, pattern='*_ANATOMY' )
print( sprintf('anatomy folders:' ) )
anatomyFolders

epiSessions <- c( epiSessions_orig, epiSessions_extra ) 

runCodeFlag <- 1

Sys.setenv(OMP_NUM_THREADS='8')
Sys.getenv('OMP_NUM_THREADS')

### align anatomyFolders to epiSession 
epiSessionSubj <- rep('aaa',length(epiSessions) )
for ( nSession in 1:length( epiSessions ) ) { # 
  currentEpiSession <- epiSessions[ nSession ]
  currentEpiSessionSplit <- strsplit( currentEpiSession, '_' )[[1]][1]
  epiSessionSubj[ nSession ] <- currentEpiSessionSplit
}

anatSubj <- rep('aaa',length(anatomyFolders) )
for ( nAnatomy in 1:length( anatomyFolders ) ) { # 
  currentAnatomy <- anatomyFolders[ nAnatomy ]
  currentAnatomySplit <- strsplit( currentAnatomy, '_' )[[1]][1]
  anatSubj[ nAnatomy ] <- currentAnatomySplit
}

anatIndex <- rep( 999, length( epiSessionSubj ) )
for( nSession in 1:length( epiSessionSubj ) ) { # nSession <- 1
  epiTemp <- epiSessionSubj[ nSession ]
  anatIndex[ nSession ] <- which( anatSubj == epiTemp )
}

#anatomyPerSession <- rep( anatomyFolders, array( table( epiSessionSubj ) ) )
anatomyPerSession <- anatomyFolders[ anatIndex ]

### check the alignment:
data.frame( anatomyPerSession, epiSessions )

### experimental conditions ###

# orig
Ses01 <- c('bars','bars','eyeBar3','eyeBar3','bars','eyeBar3','bars','eyeBar3','bars')
Ses02 <- c('eyeBar2','eyeBar2','kastnerBaseline','eyeBar2','eyeBar2','kastnerBaseline','eyeBar2','kastnerBaseline','kastnerBaseline')
Ses03 <- c('bars','bars','bars','bars','eyeSimple','eyeSimple','eyeSimple','eyeSimple','eyeSimple')
Ses04 <- c('kastnerFinerCw','kastnerFinerCCw','kastnerFinerCw','kastnerFinerCCw','vestLocalizer','vestLocalizer','vestLocalizer','kastnerFinerCw','kastnerFinerCCw')
Ses05 <- c('eyeBar3','eyeBar3','kastnerBaseline','kastnerBaseline','eyeBar3','kastnerBaseline','eyeBar3','kastnerBaseline','kastnerBaseline')
Ses06 <- c('bars','bars','bars','bars','eyeSimple','eyeSimple','eyeSimple','eyeSimple')
Ses07 <- c('kastnerFinerCw','kastnerFinerCCw','kastnerFinerCw','kastnerFinerCCw','vestLocalizer','vestLocalizer','vestLocalizer','kastnerFinerCw','kastnerFinerCCw')
Ses08 <- c('bars','bars','bars','kastnerFinerCw','kastnerFinerCCw','kastnerFinerCw','kastnerFinerCCw','bars')
Ses09 <- c('kastnerBaseline','kastnerBaseline','kastnerBaseline','kastnerBaseline','kastnerBaseline','kastnerFinerCw','kastnerFinerCCw','kastnerBaseline','kastnerBaseline')
Ses10 <- c('bars','bars','bars','kastnerFinerCw','kastnerFinerCCw','kastnerFinerCw','kastnerFinerCCw','bars')
Ses11 <- c('kastnerBaseline','kastnerBaseline','kastnerBaseline','kastnerBaseline','kastnerBaseline','kastnerFinerCw','kastnerFinerCCw','kastnerBaseline')
Ses12 <- c('bars','bars','eyeBar3','eyeBar3','bars','bars','eyeBar3','eyeBar3','bars')
Ses13 <- c('eyeBar3','eyeBar3','kastnerBaseline','kastnerBaseline','eyeBar3','eyeBar3','kastnerBaseline','kastnerBaseline')
Ses14 <- c('eyeBar2','eyeBar2','eyeBar2','eyeBar2','eyeBar2','eyeBar2')
Ses15 <- c('bars','bars','bars','bars','eyeSimple','eyeSimple','eyeSimple','eyeSimple')
Ses16 <- c('kastnerFinerCw','kastnerFinerCCw','vestLocalizer','kastnerFinerCw','kastnerFinerCCw','vestLocalizer','kastnerFinerCw','kastnerFinerCCw','vestLocalizer')
Ses17 <- c('bars','bars','kastnerFinerCw','kastnerFinerCCw','bars','kastnerFinerCw','bars','kastnerFinerCCw')
Ses18 <- c('kastnerBaseline','kastnerBaseline','kastnerBaseline','kastnerBaseline','kastnerFinerCw','kastnerFinerCCw','kastnerBaseline','kastnerBaseline','kastnerBaseline')
Ses19 <- c('bars','bars','bars','eyeSimple','eyeSimple','eyeSimple','eyeSimple','bars')
Ses20 <- c('eyeBar3','eyeBar3','kastnerBaseline','eyeBar3','kastnerBaseline','eyeBar3','kastnerBaseline','kastnerBaseline')
Ses21 <- c('kastnerFinerCw','kastnerFinerCCw','kastnerFinerCw','kastnerFinerCCw','vestLocalizer','vestLocalizer','vestLocalizer') # first epi, Kastner Finer cw, add the extra 6 TRs from the next run
Ses22 <- c('bars','bars','bars','kastnerFinerCw','kastnerFinerCCw','kastnerFinerCw','kastnerFinerCCw','bars')
Ses23 <- c('kastnerBaseline','kastnerBaseline','kastnerBaseline','kastnerBaseline','kastnerBaseline','kastnerFinerCw','kastnerFinerCCw','kastnerBaseline','kastnerBaseline')
Ses24 <- c('eyeBar2','eyeBar2','eyeBar2','eyeBar2','eyeBar2','kastnerBaseline','kastnerBaseline','kastnerBaseline','kastnerBaseline')
Ses25 <- c('bars','bars','eyeBar3','eyeBar3','bars','bars','eyeBar3','eyeBar3','bars')
Ses26 <- c('eyeBar2','eyeBar2','eyeBar2','eyeBar2','eyeBar2','kastnerFinerCw','kastnerFinerCw') # first Kastner Finer cw, pad with the second
Ses27 <- c('kastnerFinerCw','kastnerFinerCCw','kastnerFinerCw','kastnerFinerCCw','vestLocalizer','vestLocalizer','vestLocalizer') 
Ses28 <- c('eyeBar3','eyeBar3','kastnerBaseline','eyeBar3','eyeBar3','kastnerBaseline','kastnerBaseline','kastnerBaseline') 
Ses29 <- c('bars','bars','bars','bars','eyeSimple','eyeSimple','eyeSimple','eyeSimple','eyeSimple','eyeSimple','bars') 
Ses30 <- c('bars','bars','eyeBar3','eyeBar3','bars','bars','eyeBar3','eyeBar3','bars','kastnerBaseline','kastnerBaseline','kastnerBaseline','kastnerBaseline') 
Ses31 <- c('eyeBar2','kastnerBaseline','eyeBar2','kastnerBaseline','eyeBar2','eyeBar2','eyeBar2','eyeBar2','eyeBar2','kastnerBaseline','kastnerBaseline') # maybe remove last kastnerBaseline ran, aborted due to heating
Ses32 <- c('eyeSimple','eyeSimple','eyeSimple','bars','bars','eyeSimple','bars','eyeSimple','bars')
Ses33 <- c('kastnerFinerCw','kastnerFinerCCw','kastnerFinerCw','kastnerFinerCCw','vestLocalizer','vestLocalizer','kastnerFinerCw','kastnerFinerCCw','vestLocalizer','vestLocalizer')
Ses34 <- c('bars','bars','eyeSimple','eyeSimple','bars','bars','eyeSimple','eyeSimple')
Ses35 <- c('kastnerFinerCw','kastnerFinerCCw','vestLocalizer','kastnerFinerCw','kastnerFinerCCw','vestLocalizer','kastnerFinerCw','kastnerFinerCCw','vestLocalizer')

# extra
Ses36 <- c('eyeBar2','eyeBar2','eyeBar3','eyeBar3','eyeBar2','eyeBar2','eyeBar3','eyeBar3','eyeBar2')
Ses37 <- c('eyeBar2','eyeBar2','eyeBar3','eyeBar3','eyeBar2','eyeBar2','eyeBar3','eyeBar3','eyeBar2')
Ses38 <- c('eyeBar2','eyeBar2','eyeBar2','eyeBar2','eyeBar2')
Ses39 <- c('bars','bars','eyeBar3','eyeBar3','bars','eyeBar3','bars','eyeBar3','bars')
Ses40 <- c('eyeBar2','eyeBar2','eyeBar2','eyeBar2','eyeBar2')
Ses41 <- c('bars','bars','eyeBar3','eyeBar3','bars','bars','eyeBar3','eyeBar3','bars')
Ses42 <- c('eyeBar2','eyeBar2','eyeBar2','eyeBar2','eyeBar3','eyeBar3','eyeBar3','eyeBar3','eyeBar2')
Ses43 <- c('bars','bars','eyeBar3','eyeBar3','bars','eyeBar3','bars','eyeBar3','bars','eyeBar3')
Ses44 <- c('eyeBar2','eyeBar2','eyeBar2','eyeBar2','eyeBar2')

# prepare full dataset to check conditions
dfOut <- c()
for ( nSession in 1:length( epiSessions ) ) { # nSession <- 1
  
  if ( nSession < 36 ) { mainDir <- sprintf('%s/%s', epiDir, epiSessions[ nSession ] ); epiDirTemp <- epiDir }
  if ( nSession >= 36 ) { mainDir <- sprintf('%s/%s', epiDir_extraData, epiSessions[ nSession ] ); epiDirTemp <- epiDir_extraData }
  anatomyDirLoop <- sprintf('%s/%s', anatomyDir, anatomyPerSession[ nSession ] )
  
  setwd( mainDir )
  getwd()
  print( getwd() )
  print( sprintf('session: %d', nSession ) )
  
  # here you can use 'data_tsCorr_denoised_distCorr_motCorr' or 'epiOnAnat' (epi interpolated on anatomy space)
  epiFiles <- list.files( path='epiOnAnat', recursive = FALSE, full.names = FALSE )
  lengthEpiFiles <- length( epiFiles )
  
  if ( nSession<10 ) { evalText <- sprintf( 'Ses0%d', nSession ) }
  if ( nSession>=10 ) { evalText <- sprintf( 'Ses%d', nSession ) }
  sesTemp <- eval( parse( text=evalText ) )
  sessionTemp <- rep( epiSessions[ nSession ], lengthEpiFiles )
  anatTemp <- rep( anatomyPerSession[ nSession ], lengthEpiFiles )
  lengthEpiFilesTemp <- rep( lengthEpiFiles, lengthEpiFiles )
  fileIdxTemp <- seq(1,lengthEpiFiles)
  epiDirArray <- rep( epiDirTemp, lengthEpiFiles )
  # check match between n epi files and n of conditions expected
  if ( length( sesTemp ) != lengthEpiFiles ) { print( sprintf('different epi n of files, check!') ); break }
  
  dfTemp <- data.frame( anatTemp, sesTemp, sessionTemp, lengthEpiFilesTemp, fileIdxTemp, epiDirArray )
  dfOut <- rbind( dfOut, dfTemp )
}

# conditions X participant table:
table( dfOut$sesTemp, dfOut$anatTemp )

# conditions X participant, comupte average time series; - complete the missing time points in datasets, control previous code
expConditions <- unique( dfOut$sesTemp )
for ( nCondition in 1:length( expConditions ) ) { # nCondition <- 1 # length( expConditions )
  condDfTemp <- subset( dfOut, sesTemp==expConditions[ nCondition ] )
  
  #### create average condition folder folder ####
  setwd(mainFolder)
  outputFolder <- sprintf( 'expConditions/%s', expConditions[ nCondition ] )
  
  sprintf('###########')
  sprintf('###########')
  print( outputFolder )
  sprintf('###########')
  sprintf('###########')
  if ( dir.exists( outputFolder ) ) {
    instr <- sprintf( 'rm -R %s',  outputFolder )
    print( instr )
    if (runCodeFlag==1) { system( instr ) }  
  }  
  if ( !dir.exists( outputFolder ) ) {
    instr <- sprintf( 'mkdir %s',  outputFolder )
    print( instr )
    if (runCodeFlag==1) { system( instr ) }  
  }         
  
  participantsTemp <- unique( condDfTemp$anatTemp )
  for ( nPart in 1:length( participantsTemp ) ) { # nPart <- 18
    partDfTemp <- subset( condDfTemp, anatTemp==participantsTemp[ nPart ] )
    participantID <- strsplit( partDfTemp[1,1], '_' )[[1]][1] 
    instrBeginning <- sprintf('3dMean -prefix %s/%s/%s_%s.nii.gz', mainFolder, outputFolder, participantID, expConditions[ nCondition ]  )
    
    sprintf('###########')
    sprintf('###########')
    print( participantID )
    sprintf('###########')
    sprintf('###########')
    
    
    # build instruction
    for ( nReps in 1:dim( partDfTemp )[1] ) { # nReps <- 1
      #setwd( epiDir )
      setwd( partDfTemp[nReps,6] )
      setwd( partDfTemp[nReps,3] )
      print( getwd() )
      
      epiFiles <- list.files( path='epiOnAnat', recursive = FALSE, full.names = FALSE  )
      instrBeginning <- paste( instrBeginning,
                               sprintf( '%s/%s/epiOnAnat/%s', partDfTemp[nReps,6], partDfTemp[nReps,3], epiFiles[ partDfTemp[nReps,5] ]  ) )
    }
    
    print( instrBeginning )
    if (runCodeFlag==1) { system( instrBeginning ) }  
    
  }
}


