rm( list=ls() ); gc();

#need to load module R, afni, and activate conda environment beforehand

#### Move to the desired folder ####
sprintf('Move to the desired folder:')
mainDir <- '/scratch/af4887/Proj_Gaia_David'
inputDir <- '/scratch/af4887/Proj_Gaia_David/EPI_timeSlicedCorrected'
outputDirSafe <- '/scratch/af4887/Proj_Gaia_David/EPI_timeSlicedCorrected_Safe'
outputDirThreat <- '/scratch/af4887/Proj_Gaia_David/EPI_timeSlicedCorrected_Threat'
stimTimingFolder <- '/scratch/af4887/Proj_Gaia_David/wetransfer_stimTiming/threat_safe'
stimTimingRunConditions <- '/scratch/af4887/Proj_Gaia_David/wetransfer_stimTiming/run_threat_safe.csv'

setwd( mainDir )
sprintf( 'Current folder: %s', print( getwd() ) )
singleSubjectFolders <- dir( inputDir )
runCodeFlag <- 1

# clean up output folder Safe, if it exists
if ( dir.exists( outputDirSafe ) ) {
  instr <- sprintf('rm -R %s', outputDirSafe )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
}
instr <- sprintf('mkdir %s', outputDirSafe )
print( instr )
if (runCodeFlag==1) { system( instr ) }

# clean up output folder Threat, if it exists
if ( dir.exists( outputDirThreat ) ) {
  instr <- sprintf('rm -R %s', outputDirThreat )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
}
instr <- sprintf('mkdir %s', outputDirThreat )
print( instr )
if (runCodeFlag==1) { system( instr ) }

# load stimTimingRunConditions file
stimTimingRunConditionsDF <- read.csv( stimTimingRunConditions, header=TRUE )
table( stimTimingRunConditionsDF$sub )
allSubjs_stimTimingRunConditionsDF <- unique( stimTimingRunConditionsDF$sub ) 

for ( nSubj in 1 : length( singleSubjectFolders )  ) {# nSubj <- 1 

  # create individual participant folder in thread and safe folders  
  instr <- sprintf('mkdir %s/%s', outputDirSafe, singleSubjectFolders[ nSubj ] )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  instr <- sprintf('mkdir %s/%s', outputDirThreat, singleSubjectFolders[ nSubj ] )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }

  # create individual participant stim timing folders in thread and safe folders  
  instr <- sprintf('mkdir %s/%s/stimTiming', outputDirSafe, singleSubjectFolders[ nSubj ] )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  instr <- sprintf('mkdir %s/%s/stimTiming', outputDirThreat, singleSubjectFolders[ nSubj ] )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  
  subjIdxFromRunConditionDF <- allSubjs_stimTimingRunConditionsDF[ nSubj ]

  sprintf( '' )
  sprintf( '' )
  sprintf( '#######################################' )
  sprintf( '#######################################' )
  sprintf('Move to the individual participant folder:')
  setwd( inputDir )
  setwd( singleSubjectFolders[ nSubj ] )  
  sprintf( 'Current folder: %s', print( getwd() ) )
  sprintf( 'subj idx from run conditions DF: %d', subjIdxFromRunConditionDF  )
  sprintf( '#######################################' )
  sprintf( '#######################################' )
  sprintf( '' )
  sprintf( '' )
  nifti_files_toProcess <- dir( pattern='*nii.gz' )

  subjDFDatasetLoop <- subset( stimTimingRunConditionsDF, sub==subjIdxFromRunConditionDF )
  
  ## copy epi data in the correct folder
  for ( nRun in 1:dim( subjDFDatasetLoop )[1] ) { # loop along the runs and copy the specific run on the specific output folder
    if ( subjDFDatasetLoop[nRun,2]=='safe' ) {
      instr <- sprintf('cp %s %s/%s/%s', nifti_files_toProcess[nRun], outputDirSafe, singleSubjectFolders[ nSubj ], nifti_files_toProcess[nRun]  )      
      print( instr )
      if (runCodeFlag==1) { system( instr ) }
    }
    if ( subjDFDatasetLoop[nRun,2]=='threat' ) {
      instr <- sprintf('cp %s %s/%s/%s', nifti_files_toProcess[nRun], outputDirThreat, singleSubjectFolders[ nSubj ], nifti_files_toProcess[nRun]  )      
      print( instr )
      if (runCodeFlag==1) { system( instr ) }
    }
  }
  
  ## copy stim time data in the correct folder (safe or threat)
  setwd( stimTimingFolder )
  subjectStimTime_safe <- dir( pattern=sprintf( '%s_%s*', singleSubjectFolders[ nSubj ], 'safe' ) )
  subjectStimTime_threat <- dir( pattern=sprintf( '%s_%s*', singleSubjectFolders[ nSubj ], 'threat' ) )
  for ( nFile in 1:length( subjectStimTime_safe ) ) { # nFile <- 1
    instr <- sprintf('cp %s %s/%s/stimTiming', subjectStimTime_safe[nFile], outputDirSafe, singleSubjectFolders[ nSubj ] )
    print( instr )
    if (runCodeFlag==1) { system( instr ) }
  }
  for ( nFile in 1:length( subjectStimTime_threat ) ) { # nFile <- 1
    instr <- sprintf('cp %s %s/%s/stimTiming', subjectStimTime_threat[nFile], outputDirThreat, singleSubjectFolders[ nSubj ] )
    print( instr )
    if (runCodeFlag==1) { system( instr ) }
  }
  
}

setwd( mainDir )
sprintf( 'done, current folder: %s', print( getwd() ) )

