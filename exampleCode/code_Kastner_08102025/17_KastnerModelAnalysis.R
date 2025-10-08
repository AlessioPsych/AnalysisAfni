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
kastnerBaselineDir <- sprintf( '%s/expConditions/kastnerBaseline', mainFolder )
anatomyDir <- sprintf( '%s/anatomies_KastnerClassic_Freesurfer', mainFolder )
modelsDir <- sprintf('%s/modelsOutput', mainFolder )

setwd( kastnerBaselineDir )
print( sprintf('current folder: %s', getwd() ) )
kastnerBaselineParticipants <- dir( kastnerBaselineDir )
print( sprintf('kastner baseline participants:' ) )
kastnerBaselineParticipants

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
anatomyFolders_matching <- rep(0,length(kastnerBaselineParticipants))
for ( nPart in 1:length( kastnerBaselineParticipants ) ) {
  partId <- strsplit( kastnerBaselineParticipants[ nPart ], '_' )[[1]][1]
  anatIdx <- which( anatomyFolders_id == partId )
  print( anatIdx )
  anatomyFolders_matching[nPart] <- anatIdx
}

# check the alignment
selectedAnatFolders <- anatomyFolders[ anatomyFolders_matching ]
data.frame( kastnerBaselineParticipants, selectedAnatFolders )

runCodeFlag <- 1

Sys.setenv(OMP_NUM_THREADS='8')
Sys.getenv('OMP_NUM_THREADS')

#### create outputFolder folder, kastnerClassic ####
dirToCheck <- sprintf('%s/kastnerBaseline', modelsDir )
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

for ( nPart in 1:length(kastnerBaselineParticipants) ) { # nPart <- 1 # length(kastnerBaselineParticipants)
  
  epiKastnerBaselineParticipantTemp <- sprintf('%s/%s', kastnerBaselineDir, kastnerBaselineParticipants[nPart] )
  anatomyTemp <- sprintf('%s/%s/FreeSeg_result/SUMA', anatomyDir, selectedAnatFolders[nPart] )
  print( epiKastnerBaselineParticipantTemp )
  print( anatomyTemp )
  participantName <- strsplit( selectedAnatFolders[nPart], '_' )[[1]][1]
  #participantName <- 'delMe'
  
  #### create outputFolder folder, kastnerClassic ####
  dirToCheck <- sprintf('%s/kastnerBaseline/%s', modelsDir, participantName )
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

  # lh ribbon  
  instr <- paste('3dresample',
                 sprintf( '-input %s/lh.ribbon.nii', anatomyTemp ),
                 sprintf( '-master %s', epiKastnerBaselineParticipantTemp ),  
                 sprintf( '-prefix %s/%s/%s/%s.zzz.lh.ribbon.nii', modelsDir, 'kastnerBaseline', participantName, selectedAnatFolders[nPart] ), 
                 '-rmode NN' )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }  
  
  # rh ribbon  
  instr <- paste('3dresample',
                 sprintf( '-input %s/rh.ribbon.nii', anatomyTemp ),
                 sprintf( '-master %s', epiKastnerBaselineParticipantTemp ),  
                 sprintf( '-prefix %s/%s/%s/%s.zzz.rh.ribbon.nii', modelsDir, 'kastnerBaseline', participantName, selectedAnatFolders[nPart] ), 
                 '-rmode NN' )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }  
  
  # combine ribbons
  instr <- paste('3dcalc', 
                 sprintf( '-a %s/%s/%s/%s.zzz.lh.ribbon.nii', modelsDir, 'kastnerBaseline', participantName, selectedAnatFolders[nPart] ), 
                 sprintf( '-b %s/%s/%s/%s.zzz.rh.ribbon.nii', modelsDir, 'kastnerBaseline', participantName, selectedAnatFolders[nPart] ), 
                 sprintf('-expr \u0027a+b\u0027'),
                 sprintf( '-prefix %s/%s/%s/%s.zzz.grayMatter.nii', modelsDir, 'kastnerBaseline', participantName, selectedAnatFolders[nPart] ) ) 
  print( instr )
  if (runCodeFlag==1) { system( instr ) } 
  
  # copy kastner baseline data
  instr <- paste('cp',
                 sprintf( '%s', epiKastnerBaselineParticipantTemp ),  
                 sprintf( '%s/%s/%s/%s', modelsDir, 'kastnerBaseline', participantName, kastnerBaselineParticipants[nPart] )
  )
  print( instr )
  if (runCodeFlag==1) { system( instr ) } 
  
  #blur kastner baseline data
  instr <- paste('3dmerge',
                 '-1blur_fwhm 1.5',
                 '-doall',
                 sprintf('-prefix %s_Kastner_Baseline_blur.nii.gz', participantName ),
                 sprintf( '%s', kastnerBaselineParticipants[nPart]  ) )
  print( instr )
  if (runCodeFlag==1) { system( instr ) } 
  
  #mean kastner baseline data
  instr <- paste('3dTstat',
                 sprintf( '-mean -prefix %s_Kastner_Baseline_mean.nii', participantName ),
                 sprintf('%s_Kastner_Baseline_blur.nii.gz', participantName ) )
  print( instr )
  if (runCodeFlag==1) { system( instr ) } 
  
  # scale kastner baseline data
  instr <- paste('3dcalc',
                 sprintf('-a %s_Kastner_Baseline_blur.nii.gz', participantName ),
                 sprintf( '-b %s_Kastner_Baseline_mean.nii', participantName ),
                 sprintf( '-expr \u0027 min(200, a/b*100)*step(a)*step(b) \u0027' ),
                 sprintf( '-prefix %s_zzz_Kastner_Baseline_scaled.nii', participantName )
  )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  
  # detrend kastner baseline
  instr <- paste('3dDetrend',
                 sprintf('-polort 3'),
                 sprintf( '-prefix %s_zzz_Kastner_Baseline_detrended.nii', participantName ),
                 sprintf( '%s_zzz_Kastner_Baseline_scaled.nii', participantName )
  )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  
  #kastner model prediction, baseline model
  instr <- paste('Rscript',
                 sprintf('%s/code/models/generatePredictionsKastner_model_PSOCK.R', mainFolder ),
                 sprintf('%s_zzz_Kastner_Baseline_detrended.nii', participantName ),
                 sprintf('%s_kastner_model_baseline_model/', participantName ),
                 sprintf('0 baselineModel') )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }

  #kastner model prediction, kinematic baseline model
  instr <- paste('Rscript',
                 sprintf('%s/code/models/generatePredictionsKastner_model_PSOCK.R', mainFolder ),
                 sprintf('%s_zzz_Kastner_Baseline_detrended.nii', participantName ),
                 sprintf('%s_kastner_model_kinematic_baselineModel/', participantName ),
                 sprintf('0 kinematic_baselineModel') )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }

  #kastner model prediction, kinematic model
  instr <- paste('Rscript',
                 sprintf('%s/code/models/generatePredictionsKastner_model_PSOCK.R', mainFolder ),
                 sprintf('%s_zzz_Kastner_Baseline_detrended.nii', participantName ),
                 sprintf('%s_kastner_kinematicModel/', participantName ),
                 sprintf('0 kinematicModel') )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  
  # get mask
  instr <- sprintf('3dmask_tool -input %s_ANATOMY.zzz.grayMatter.nii -prefix %s_mask_auto.nii.gz -dilate_input 1', participantName, participantName )  
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  
  #kastner model baseline model
  instr <- paste('Rscript',
                 sprintf('%s/code/models/fitEyeKinematics_PSOCK.R', mainFolder), 
                 sprintf('%s_zzz_Kastner_Baseline_detrended.nii', participantName ), 
                 sprintf('%s_mask_auto.nii.gz', participantName ), 
                 sprintf('%s_kastner_model_baseline_model/', participantName ),
                 sprintf('%s_kastner_model_baseline', participantName ), 
                 sprintf('baselineModel') )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }

  #kastner model kinematic baseline model
  instr <- paste('Rscript',
                 sprintf('%s/code/models/fitEyeKinematics_PSOCK.R', mainFolder), 
                 sprintf('%s_zzz_Kastner_Baseline_detrended.nii', participantName ), 
                 sprintf('%s_mask_auto.nii.gz', participantName ), 
                 sprintf('%s_kastner_model_kinematic_baselineModel/', participantName ),
                 sprintf('%s_kastner_kinematic_baseline', participantName ), 
                 sprintf('kinematic_baselineModel') )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  
  #kastner model kinematic model
  instr <- paste('Rscript',
                 sprintf('%s/code/models/fitEyeKinematics_PSOCK.R', mainFolder), 
                 sprintf('%s_zzz_Kastner_Baseline_detrended.nii', participantName ), 
                 sprintf('%s_mask_auto.nii.gz', participantName ), 
                 sprintf('%s_kastner_kinematicModel/', participantName ),
                 sprintf('%s_kastner_kinematicModel', participantName ), 
                 sprintf('kinematicModel') )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  

  # clean up
  setwd( dirToCheck )
  print( getwd() )
  system('rm *zzz*')
  system( sprintf('rm %s', kastnerBaselineParticipants[nPart] ) )
  system( sprintf('rm *_blur.nii.gz' ) )

}

