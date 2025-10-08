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
mainFolder <- '/media/alessiofracasso/DATADRIVE1/KastnerData/staging_area_Kastner'
cwDir <- sprintf( '%s/expConditions/kastnerFinerCw', mainFolder )
ccwDir <- sprintf( '%s/expConditions/kastnerFinerCCw', mainFolder )
anatomyDir <- sprintf( '%s/anatomies_KastnerClassic_Freesurfer', mainFolder )
modelsPrfDir <- sprintf('%s/modelsOutput/pRF', mainFolder )
modelsDir <- sprintf('%s/modelsOutput', mainFolder )
codeFolder <- sprintf('%s/code/models', mainFolder )

setwd( cwDir )
print( sprintf('current folder: %s', getwd() ) )
cwParticipants <- dir( cwDir )
print( sprintf('cw participants:' ) )
cwParticipants

setwd( ccwDir )
print( sprintf('current folder: %s', getwd() ) )
ccwParticipants <- dir( ccwDir )
print( sprintf('ccw participants:' ) )
ccwParticipants

setwd( modelsPrfDir )
print( sprintf('current folder: %s', getwd() ) )
pRFParticipants <- dir( modelsPrfDir)
print( sprintf('pRF participants:' ) )
pRFParticipants

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
anatomyFolders_matching <- rep(0,length(cwParticipants))
for ( nPart in 1:length( cwParticipants ) ) {
  cwPartId <- strsplit( cwParticipants[ nPart ], '_' )[[1]][1]
  anatIdx <- which( anatomyFolders_id == cwPartId )
  print( anatIdx )
  anatomyFolders_matching[nPart] <- anatIdx
}

# look for corresponding pRF folders
pRFFolders_id <- rep('a',length( pRFParticipants ))
for (nPrf in 1:length(  pRFParticipants  )) {
  pRFFolders_id[nPrf] <- strsplit( pRFParticipants, '_' )[[nPrf]][1]
}
pRFFolders_matching <- rep(0,length(cwParticipants))
for ( nPart in 1:length( cwParticipants ) ) {
  cwPartId <- strsplit( cwParticipants[ nPart ], '_' )[[1]][1]
  pRFIdx <- which( pRFFolders_id == cwPartId )
  print( pRFIdx )
  pRFFolders_matching[nPart] <- pRFIdx
}

# check the alignment
selectedAnatFolders <- anatomyFolders[ anatomyFolders_matching ]
selectedpRFFolders <- pRFParticipants[ pRFFolders_matching ]
data.frame( cwParticipants, ccwParticipants, selectedAnatFolders, selectedpRFFolders )

runCodeFlag <- 1

#Sys.setenv(OMP_NUM_THREADS='8')
#Sys.getenv('OMP_NUM_THREADS')

#### create outputFolder folder, kastnerClassic_contrastBased ####
dirToCheck <- sprintf('%s/kastnerClassic_contrastBased', modelsDir )
print( dirToCheck )
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

for ( nPart in 1:length(cwParticipants) ) { # length(cwParticipants) # nPart <- 1

  epiCwTemp <- sprintf('%s/%s', cwDir, cwParticipants[nPart] )
  epiCCwTemp <- sprintf('%s/%s', ccwDir, ccwParticipants[nPart] )
  dspRFTemp <- sprintf('%s/%s/%s/%s_model_wAverage_05.nii.gz', modelsDir, 'kastnerClassic', selectedpRFFolders[nPart], selectedpRFFolders[nPart] )
  anatomyTemp <- sprintf('%s/%s/FreeSeg_result/SUMA', anatomyDir, selectedAnatFolders[nPart] )
  pRFTemp <- sprintf('%s/%s/%s_pRF_params.nii.gz', modelsPrfDir, selectedpRFFolders[nPart], selectedpRFFolders[nPart] )
  print( epiCwTemp )
  print( epiCCwTemp )
  print( dspRFTemp )
  print( anatomyTemp )
  print( pRFTemp )
  participantName <- strsplit( selectedAnatFolders[nPart], '_' )[[1]][1]
  #participantName <- 'delMe'
  
  #### create outputFolder folder, kastnerClassic ####
  dirToCheck <- sprintf('%s/kastnerClassic_contrastBased/%s', modelsDir, participantName )
  print( dirToCheck )
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
  
  setwd( dirToCheck )
  print( getwd() )  
  
  # create stimulus cw
  if ( !file.exists( 'cw_completeStimuliKastnerOrig_borderOnlyAfterSaccade_shifted_final.Rdata' ) ) {
    instr <- sprintf( 'Rscript %s/contrastBasedExpectationFromPRF_buildStimuli.R cw', codeFolder )
    print( instr )
    if (runCodeFlag==1) { system( instr ) }
  }
  
  # create stimulus ccw
  if ( !file.exists( 'ccw_completeStimuliKastnerOrig_borderOnlyAfterSaccade_shifted_final.Rdata' ) ) {
    instr <- sprintf( 'Rscript %s/contrastBasedExpectationFromPRF_buildStimuli.R ccw', codeFolder )
    print( instr )
    if (runCodeFlag==1) { system( instr ) }
  }
  
  # fit stimulus cw
  instr <- sprintf( 'Rscript %s/contrastBasedExpectationFromPRF_fit_parallel02.R %s %s %s %s %s %s %s', 
                    codeFolder, 
                    epiCwTemp,
                    pRFTemp,
                    'cw_completeStimuliKastnerOrig_borderOnlyAfterSaccade_shifted_final.Rdata',
                    'cw',
                    dirToCheck,
                    participantName,
                    dspRFTemp )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  
  # fit stimulus ccw
  instr <- sprintf( 'Rscript %s/contrastBasedExpectationFromPRF_fit_parallel02.R %s %s %s %s %s %s %s', 
                    codeFolder, 
                    epiCCwTemp,
                    pRFTemp,
                    'ccw_completeStimuliKastnerOrig_borderOnlyAfterSaccade_shifted_final.Rdata',
                    'ccw',
                    dirToCheck,
                    participantName,
                    dspRFTemp )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  
}

