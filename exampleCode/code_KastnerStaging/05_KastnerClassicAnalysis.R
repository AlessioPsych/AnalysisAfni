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
cwDir <- sprintf( '%s/expConditions/kastnerFinerCw', mainFolder )
ccwDir <- sprintf( '%s/expConditions/kastnerFinerCCw', mainFolder )
anatomyDir <- sprintf( '%s/anatomies_KastnerClassic_Freesurfer', mainFolder )
modelsDir <- sprintf('%s/modelsOutput', mainFolder )

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

# check the alignment
selectedAnatFolders <- anatomyFolders[ anatomyFolders_matching ]
data.frame( cwParticipants, ccwParticipants, selectedAnatFolders )

runCodeFlag <- 1

Sys.setenv(OMP_NUM_THREADS='8')
Sys.getenv('OMP_NUM_THREADS')

#### create outputFolder folder, kastnerClassic ####
dirToCheck <- sprintf('%s/kastnerClassic', modelsDir )
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

for ( nPart in 1:length(cwParticipants) ) { #length(cwParticipants) # nPart <- 1
  
  epiCwTemp <- sprintf('%s/%s', cwDir, cwParticipants[nPart] )
  epiCCwTemp <- sprintf('%s/%s', ccwDir, ccwParticipants[nPart] )
  anatomyTemp <- sprintf('%s/%s/FreeSeg_result/SUMA', anatomyDir, selectedAnatFolders[nPart] )
  print( epiCwTemp )
  print( epiCCwTemp )
  print( anatomyTemp )
  participantName <- strsplit( selectedAnatFolders[nPart], '_' )[[1]][1]
  #participantName <- 'delMe'
  
  #### create outputFolder folder, kastnerClassic ####
  dirToCheck <- sprintf('%s/kastnerClassic/%s', modelsDir, participantName )
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
                 sprintf( '-master %s', epiCwTemp ),  
                 sprintf( '-prefix %s/%s/%s/%s.zzz.lh.ribbon.nii', modelsDir, 'kastnerClassic', participantName, selectedAnatFolders[nPart] ), 
                 '-rmode NN' )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }  
  
  # rh ribbon  
  instr <- paste('3dresample',
                 sprintf( '-input %s/rh.ribbon.nii', anatomyTemp ),
                 sprintf( '-master %s', epiCwTemp ),  
                 sprintf( '-prefix %s/%s/%s/%s.zzz.rh.ribbon.nii', modelsDir, 'kastnerClassic', participantName, selectedAnatFolders[nPart] ), 
                 '-rmode NN' )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }  
  
  # combine ribbons
  instr <- paste('3dcalc', 
                 sprintf( '-a %s/%s/%s/%s.zzz.lh.ribbon.nii', modelsDir, 'kastnerClassic', participantName, selectedAnatFolders[nPart] ), 
                 sprintf( '-b %s/%s/%s/%s.zzz.rh.ribbon.nii', modelsDir, 'kastnerClassic', participantName, selectedAnatFolders[nPart] ), 
                 sprintf('-expr \u0027a+b\u0027'),
                 sprintf( '-prefix %s/%s/%s/%s.zzz.grayMatter.nii', modelsDir, 'kastnerClassic', participantName, selectedAnatFolders[nPart] ) ) 
  print( instr )
  if (runCodeFlag==1) { system( instr ) } 
  
  # copy cw data
  instr <- paste('cp',
                 sprintf( '%s', epiCwTemp ),  
                 sprintf( '%s/%s/%s/%s', modelsDir, 'kastnerClassic', participantName, cwParticipants[nPart] )
  )
  print( instr )
  if (runCodeFlag==1) { system( instr ) } 
  
  # copy ccw data
  instr <- paste('cp',
                 sprintf( '%s', epiCCwTemp ),  
                 sprintf( '%s/%s/%s/%s', modelsDir, 'kastnerClassic', participantName, ccwParticipants[nPart] )
  )
  print( instr )
  if (runCodeFlag==1) { system( instr ) } 
  
  #blur cw
  instr <- paste('3dmerge',
                 '-1blur_fwhm 1.5',
                 '-doall',
                 sprintf('-prefix %s_CW_blur.nii.gz', participantName ),
                 sprintf( '%s', cwParticipants[nPart]  ) )
  print( instr )
  if (runCodeFlag==1) { system( instr ) } 
  
  #blur ccw
  instr <- paste('3dmerge',
                 '-1blur_fwhm 1.5',
                 '-doall',
                 sprintf('-prefix %s_CCW_blur.nii.gz', participantName ),
                 sprintf( '%s', ccwParticipants[nPart]  ) )
  print( instr )
  if (runCodeFlag==1) { system( instr ) } 
  
  #mean cw
  instr <- paste('3dTstat',
                 sprintf( '-mean -prefix %s_CW_mean.nii', participantName ),
                 sprintf('%s_CW_blur.nii.gz', participantName ) )
  print( instr )
  if (runCodeFlag==1) { system( instr ) } 
  
  #mean ccw
  instr <- paste('3dTstat',
                 sprintf( '-mean -prefix %s_CCW_mean.nii', participantName ),
                 sprintf('%s_CCW_blur.nii.gz', participantName ) )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  
  # scale cw
  instr <- paste('3dcalc',
                 sprintf('-a %s_CW_blur.nii.gz', participantName ),
                 sprintf( '-b %s_CW_mean.nii', participantName ),
                 sprintf( '-expr \u0027 min(200, a/b*100)*step(a)*step(b) \u0027' ),
                 sprintf( '-prefix %s_zzz_CW_scaled.nii', participantName )
  )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  
  # scale ccw
  instr <- paste('3dcalc',
                 sprintf('-a %s_CCW_blur.nii.gz', participantName ),
                 sprintf( '-b %s_CCW_mean.nii', participantName ),
                 sprintf( '-expr \u0027 min(200, a/b*100)*step(a)*step(b) \u0027' ),
                 sprintf( '-prefix %s_zzz_CCW_scaled.nii', participantName )
  )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  
  # detrend cw
  instr <- paste('3dDetrend',
                 sprintf('-polort 4'),
                 sprintf( '-prefix %s_zzz_CW_detrended.nii', participantName ),
                 sprintf( '%s_zzz_CW_scaled.nii', participantName )
  )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  
  # detrend ccw
  instr <- paste('3dDetrend',
                 sprintf('-polort 4'),
                 sprintf( '-prefix %s_zzz_CCW_detrended.nii', participantName ),
                 sprintf( '%s_zzz_CCW_scaled.nii', participantName )
  )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  
  #travelling wave, 05 cycles, CW
  instr <- sprintf('compTravellingWave.sh %s_zzz_CW_detrended.nii 5', participantName)
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  instr <- sprintf('mv trWave.nii.gz %s_trWave05_CW.nii.gz', participantName )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  instr <- sprintf('mv periodogram.nii.gz %s_periodogram05_CW.nii.gz', participantName )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  
  #travelling wave, 05 cycles, CCW
  instr <- sprintf('compTravellingWave.sh %s_zzz_CCW_detrended.nii 5', participantName)
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  instr <- sprintf('mv trWave.nii.gz %s_trWave05_CCW.nii.gz', participantName )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  instr <- sprintf('mv periodogram.nii.gz %s_periodogram05_CCW.nii.gz', participantName )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  
  #travelling wave, 10 cycles, CW
  instr <- sprintf('compTravellingWave.sh %s_zzz_CW_detrended.nii 10', participantName)
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  instr <- sprintf('mv trWave.nii.gz %s_trWave10_CW.nii.gz', participantName )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  instr <- sprintf('mv periodogram.nii.gz %s_periodogram10_CW.nii.gz', participantName )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  
  #travelling wave, 10 cycles, CCW
  instr <- sprintf('compTravellingWave.sh %s_zzz_CCW_detrended.nii 10', participantName)
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  instr <- sprintf('mv trWave.nii.gz %s_trWave10_CCW.nii.gz', participantName )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  instr <- sprintf('mv periodogram.nii.gz %s_periodogram10_CCW.nii.gz', participantName )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  
  # (weighted) average phase-encoded, shifted results, 05 cycles:
  CWmod <- read.AFNI( sprintf( '%s_trWave05_CW.nii.gz', participantName ) )
  CCWmod <- read.AFNI( sprintf( '%s_trWave05_CCW.nii.gz', participantName ) ) # this is the volume (CCW) where right hemi phase values are [0,pi] and left hemi value are [pi,2pi], like the modelling results, so keep this volume as the base and shift the phase values of the other volume (CW)
  CWmodVol <- CWmod$brk
  CCWmodVol <- CCWmod$brk
  outVol <- array(0, dim( CCWmod$brk ) )
  weightVolumeCW <- CWmodVol[,,,2]
  weightVolumeCCW <- CCWmodVol[,,,2]
  for (nVols in 1:dim( CCWmod$brk )[4] ) {
    if (nVols!=3) { 
      outVol[,,,nVols] <- ( CWmodVol[,,,nVols]*weightVolumeCW + CCWmodVol[,,,nVols]*weightVolumeCCW ) / ( weightVolumeCW + weightVolumeCCW ) 
    }
    if (nVols==3) {
      xCartVolCW <- cos( abs( CWmodVol[,,,nVols]-2*pi ) - pi ) # convert to cartesian assuming a radius or 1, and shifting CW values to match CCW values --> abs(x-2*pi)
      yCartVolCW <- sin( abs( CWmodVol[,,,nVols]-2*pi ) - pi ) # convert to cartesian assuming a radius or 1, and shifting CW values to match CCW values --> abs(x-2*pi)
      xCartVolCCW <- cos( CCWmodVol[,,,nVols]-pi ) # convert to cartesian assuming a radius or 1
      yCartVolCCW <- sin( CCWmodVol[,,,nVols]-pi ) # convert to cartesian assuming a radius or 1
      xVolWeighted <- ( xCartVolCW*weightVolumeCW + xCartVolCCW*weightVolumeCCW ) / ( weightVolumeCW + weightVolumeCCW )
      yVolWeighted <- ( yCartVolCW*weightVolumeCW + yCartVolCCW*weightVolumeCCW ) / ( weightVolumeCW + weightVolumeCCW )
      outVol[,,,nVols] <- atan2( yVolWeighted, xVolWeighted ) + pi
    }
  }
  write.AFNI( sprintf('%s_phase_encoded_wAverage_05.nii.gz', participantName ), brk=outVol, origin=CWmod$origin, orient=CWmod$orient, defhead=CWmod$NI_head )
  # shift output, atan2 gives the inverted output
  instr <- sprintf('compShiftPhase_model.sh %s_phase_encoded_wAverage_05.nii.gz 180 %s_phase_encoded_wAverage_05_correct.nii.gz 3', participantName, participantName )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  
  # (weighted) average phase-encoded, shifted results, 10 cycles:
  CWmod <- read.AFNI( sprintf( '%s_trWave10_CW.nii.gz', participantName ) )
  CCWmod <- read.AFNI( sprintf( '%s_trWave10_CCW.nii.gz', participantName ) ) # this is the volume (CCW) where right hemi phase values are [0,pi] and left hemi value are [pi,2pi], like the modelling results, so keep this volume as the base and shift the phase values of the other volume (CW)
  CWmodVol <- CWmod$brk
  CCWmodVol <- CCWmod$brk
  outVol <- array(0, dim( CCWmod$brk ) )
  weightVolumeCW <- CWmodVol[,,,2]
  weightVolumeCCW <- CCWmodVol[,,,2]
  for (nVols in 1:dim( CCWmod$brk )[4] ) {
    if (nVols!=3) { 
      outVol[,,,nVols] <- ( CWmodVol[,,,nVols]*weightVolumeCW + CCWmodVol[,,,nVols]*weightVolumeCCW ) / ( weightVolumeCW + weightVolumeCCW ) 
    }
    if (nVols==3) {
      xCartVolCW <- cos( abs( CWmodVol[,,,nVols]-2*pi ) - pi ) # convert to cartesian assuming a radius or 1, and shifting CW values to match CCW values --> abs(x-2*pi)
      yCartVolCW <- sin( abs( CWmodVol[,,,nVols]-2*pi ) - pi ) # convert to cartesian assuming a radius or 1, and shifting CW values to match CCW values --> abs(x-2*pi)
      xCartVolCCW <- cos( CCWmodVol[,,,nVols] ) # convert to cartesian assuming a radius or 1
      yCartVolCCW <- sin( CCWmodVol[,,,nVols] ) # convert to cartesian assuming a radius or 1
      xVolWeighted <- ( xCartVolCW*weightVolumeCW + xCartVolCCW*weightVolumeCCW ) / ( weightVolumeCW + weightVolumeCCW )
      yVolWeighted <- ( yCartVolCW*weightVolumeCW + yCartVolCCW*weightVolumeCCW ) / ( weightVolumeCW + weightVolumeCCW )
      outVol[,,,nVols] <- atan2( yVolWeighted, xVolWeighted ) + pi
    }
  }
  write.AFNI( sprintf('%s_phase_encoded_wAverage_10.nii.gz', participantName), brk=outVol, origin=CWmod$origin, orient=CWmod$orient, defhead=CWmod$NI_head )
  # shift output, atan2 gives the inverted output
  instr <- sprintf('compShiftPhase_model.sh %s_phase_encoded_wAverage_10.nii.gz 180 %s_phase_encoded_wAverage_10_correct.nii.gz 3', participantName, participantName )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  
  #kastner model prediction 2 components CW
  instr <- paste('Rscript',
                 sprintf('%s/code/models/generatePredictionsKastner_2Components_differentDirection_PSOCK.R', mainFolder ),
                 sprintf('%s_zzz_CW_detrended.nii', participantName ),
                 sprintf('%s_kastnerCW_2components_both/', participantName ),
                 sprintf('0 cw both') )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  
  #kastner model prediction 2 components CCW
  instr <- paste('Rscript',
                 sprintf('%s/code/models/generatePredictionsKastner_2Components_differentDirection_PSOCK.R', mainFolder ),
                 sprintf('%s_zzz_CCW_detrended.nii', participantName ),
                 sprintf('%s_kastnerCCW_2components_both/', participantName ),
                 sprintf('0 ccw both') )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  
  #kastner model prediction 2 components CW 10
  instr <- paste('Rscript',
                 sprintf('%s/code/models/generatePredictionsKastner_2Components_differentDirection_PSOCK.R', mainFolder ),
                 sprintf('%s_zzz_CW_detrended.nii', participantName ),
                 sprintf('%s_kastnerCW_2components_10/', participantName ),
                 sprintf('0 cw 10') )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  
  #kastner model prediction 2 components CCW 10
  instr <- paste('Rscript',
                 sprintf('%s/code/models/generatePredictionsKastner_2Components_differentDirection_PSOCK.R', mainFolder ),
                 sprintf('%s_zzz_CCW_detrended.nii', participantName ),
                 sprintf('%s_kastnerCCW_2components_10/', participantName ),
                 sprintf('0 ccw 10') )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  
  #kastner model prediction 2 components CW 5
  instr <- paste('Rscript',
                 sprintf('%s/code/models/generatePredictionsKastner_2Components_differentDirection_PSOCK.R', mainFolder ),
                 sprintf('%s_zzz_CW_detrended.nii', participantName ),
                 sprintf('%s_kastnerCW_2components_05/', participantName ),
                 sprintf('0 cw 5') )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  
  #kastner model prediction 2 components CCW 5
  instr <- paste('Rscript',
                 sprintf('%s/code/models/generatePredictionsKastner_2Components_differentDirection_PSOCK.R', mainFolder ),
                 sprintf('%s_zzz_CCW_detrended.nii', participantName ),
                 sprintf('%s_kastnerCCW_2components_05/', participantName ),
                 sprintf('0 ccw 5') )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  
  # get mask
  instr <- sprintf('3dmask_tool -input %s_ANATOMY.zzz.grayMatter.nii -prefix %s_mask_auto.nii.gz -dilate_input 1', participantName, participantName )  
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  
  #kastner model 2 components, CW both
  instr <- paste('Rscript',
                 sprintf('%s/code/models/fitEyeKinematics_PSOCK.R', mainFolder), 
                 sprintf('%s_zzz_CW_detrended.nii', participantName ), 
                 sprintf('%s_mask_auto.nii.gz', participantName ), 
                 sprintf('%s_kastnerCW_2components_both/', participantName ),
                 sprintf('%s_CW_2components_both', participantName ), 
                 sprintf('kinematicModel') )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }

  #kastner model 2 components, CCW both
  instr <- paste('Rscript',
                 sprintf('%s/code/models/fitEyeKinematics_PSOCK.R', mainFolder), 
                 sprintf('%s_zzz_CCW_detrended.nii', participantName ), 
                 sprintf('%s_mask_auto.nii.gz', participantName ), 
                 sprintf('%s_kastnerCCW_2components_both/', participantName ),
                 sprintf('%s_CCW_2components_both', participantName ), 
                 sprintf('kinematicModel') )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }

  #kastner model 2 components, CW 10
  instr <- paste('Rscript',
                 sprintf('%s/code/models/fitEyeKinematics_PSOCK.R', mainFolder), 
                 sprintf('%s_zzz_CW_detrended.nii', participantName ), 
                 sprintf('%s_mask_auto.nii.gz', participantName ), 
                 sprintf('%s_kastnerCW_2components_10/', participantName ),
                 sprintf('%s_CW_2components_10', participantName ), 
                 sprintf('kinematicModel') )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  
  #kastner model 2 components, CCW 10
  instr <- paste('Rscript',
                 sprintf('%s/code/models/fitEyeKinematics_PSOCK.R', mainFolder), 
                 sprintf('%s_zzz_CCW_detrended.nii', participantName ), 
                 sprintf('%s_mask_auto.nii.gz', participantName ), 
                 sprintf('%s_kastnerCCW_2components_10/', participantName ),
                 sprintf('%s_CCW_2components_10', participantName ), 
                 sprintf('kinematicModel') )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }

  #kastner model 2 components, CW 05
  instr <- paste('Rscript',
                 sprintf('%s/code/models/fitEyeKinematics_PSOCK.R', mainFolder), 
                 sprintf('%s_zzz_CW_detrended.nii', participantName ), 
                 sprintf('%s_mask_auto.nii.gz', participantName ), 
                 sprintf('%s_kastnerCW_2components_05/', participantName ),
                 sprintf('%s_CW_2components_05', participantName ), 
                 sprintf('kinematicModel') )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  
  #kastner model 2 components, CCW 05
  instr <- paste('Rscript',
                 sprintf('%s/code/models/fitEyeKinematics_PSOCK.R', mainFolder), 
                 sprintf('%s_zzz_CCW_detrended.nii', participantName ), 
                 sprintf('%s_mask_auto.nii.gz', participantName ), 
                 sprintf('%s_kastnerCCW_2components_05/', participantName ),
                 sprintf('%s_CCW_2components_05', participantName ), 
                 sprintf('kinematicModel') )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  
  # (weighted) average modelling results, two components:
  CWmod <- read.AFNI( sprintf('%s_CW_2components_both_params.nii.gz', participantName) )
  CCWmod <- read.AFNI( sprintf('%s_CCW_2components_both_params.nii.gz', participantName) )
  CWmodVol <- CWmod$brk
  CCWmodVol <- CCWmod$brk
  outVol <- array(0, dim( CCWmod$brk ) )
  weightVolumeCW <- CWmodVol[,,,9]
  weightVolumeCCW <- CCWmodVol[,,,9]
  for (nVols in 1:dim( CCWmod$brk )[4] ) {
    if (nVols!=6) { 
      outVol[,,,nVols] <- ( CWmodVol[,,,nVols]*weightVolumeCW + CCWmodVol[,,,nVols]*weightVolumeCCW ) / ( weightVolumeCW + weightVolumeCCW ) 
    }
    if (nVols==6) {
      xCartVolCW <- cos( CWmodVol[,,,nVols]-pi ) # convert to cartesian assuming a radius or 1
      yCartVolCW <- sin( CWmodVol[,,,nVols]-pi ) # convert to cartesian assuming a radius or 1
      xCartVolCCW <- cos( CCWmodVol[,,,nVols]-pi ) # convert to cartesian assuming a radius or 1
      yCartVolCCW <- sin( CCWmodVol[,,,nVols]-pi ) # convert to cartesian assuming a radius or 1
      xVolWeighted <- ( xCartVolCW*weightVolumeCW + xCartVolCCW*weightVolumeCCW ) / ( weightVolumeCW + weightVolumeCCW )
      yVolWeighted <- ( yCartVolCW*weightVolumeCW + yCartVolCCW*weightVolumeCCW ) / ( weightVolumeCW + weightVolumeCCW )
      outVol[,,,nVols] <- atan2( yVolWeighted, xVolWeighted ) + pi
    }
  }
  write.AFNI( sprintf( '%s_model_wAverage_both.nii.gz', participantName ), brk=outVol, origin=CWmod$origin, orient=CWmod$orient, defhead=CWmod$NI_head )
  instr <- sprintf( 'compShiftPhase_model.sh %s_model_wAverage_both.nii.gz 180 %s_model_wAverage_both_correct.nii.gz 6', participantName, participantName )
  
  # (weighted) average modelling results, 05 cycles:
  CWmod <- read.AFNI( sprintf('%s_CW_2components_05_params.nii.gz', participantName) )
  CCWmod <- read.AFNI( sprintf('%s_CCW_2components_05_params.nii.gz', participantName) )
  CWmodVol <- CWmod$brk
  CCWmodVol <- CCWmod$brk
  outVol <- array(0, dim( CCWmod$brk ) )
  weightVolumeCW <- CWmodVol[,,,9]
  weightVolumeCCW <- CCWmodVol[,,,9]
  for (nVols in 1:dim( CCWmod$brk )[4] ) {
    if (nVols!=6) { 
      outVol[,,,nVols] <- ( CWmodVol[,,,nVols]*weightVolumeCW + CCWmodVol[,,,nVols]*weightVolumeCCW ) / ( weightVolumeCW + weightVolumeCCW ) 
    }
    if (nVols==6) {
      xCartVolCW <- cos( CWmodVol[,,,nVols]-pi ) # convert to cartesian assuming a radius or 1
      yCartVolCW <- sin( CWmodVol[,,,nVols]-pi ) # convert to cartesian assuming a radius or 1
      xCartVolCCW <- cos( CCWmodVol[,,,nVols]-pi ) # convert to cartesian assuming a radius or 1
      yCartVolCCW <- sin( CCWmodVol[,,,nVols]-pi ) # convert to cartesian assuming a radius or 1
      xVolWeighted <- ( xCartVolCW*weightVolumeCW + xCartVolCCW*weightVolumeCCW ) / ( weightVolumeCW + weightVolumeCCW )
      yVolWeighted <- ( yCartVolCW*weightVolumeCW + yCartVolCCW*weightVolumeCCW ) / ( weightVolumeCW + weightVolumeCCW )
      outVol[,,,nVols] <- atan2( yVolWeighted, xVolWeighted ) + pi
    }
  }
  write.AFNI( sprintf( '%s_model_wAverage_05.nii.gz', participantName ), brk=outVol, origin=CWmod$origin, orient=CWmod$orient, defhead=CWmod$NI_head )
  instr <- sprintf('compShiftPhase_model.sh %s_model_wAverage_05.nii.gz 180 %s_model_wAverage_05_correct.nii.gz 6', participantName, participantName )
  
  # (weighted) average modelling results, 10 cycles:
  CWmod <- read.AFNI( sprintf('%s_CW_2components_10_params.nii.gz', participantName) )
  CCWmod <- read.AFNI( sprintf('%s_CCW_2components_10_params.nii.gz', participantName) )
  CWmodVol <- CWmod$brk
  CCWmodVol <- CCWmod$brk
  outVol <- array(0, dim( CCWmod$brk ) )
  weightVolumeCW <- CWmodVol[,,,9]
  weightVolumeCCW <- CCWmodVol[,,,9]
  for (nVols in 1:dim( CCWmod$brk )[4] ) {
    if (nVols!=6) { 
      outVol[,,,nVols] <- ( CWmodVol[,,,nVols]*weightVolumeCW + CCWmodVol[,,,nVols]*weightVolumeCCW ) / ( weightVolumeCW + weightVolumeCCW ) 
    }
    if (nVols==6) {
      xCartVolCW <- cos( CWmodVol[,,,nVols] - pi ) # convert to cartesian assuming a radius or 1
      yCartVolCW <- sin( CWmodVol[,,,nVols] - pi ) # convert to cartesian assuming a radius or 1
      xCartVolCCW <- cos( CCWmodVol[,,,nVols] - pi ) # convert to cartesian assuming a radius or 1
      yCartVolCCW <- sin( CCWmodVol[,,,nVols] - pi ) # convert to cartesian assuming a radius or 1
      xVolWeighted <- ( xCartVolCW*weightVolumeCW + xCartVolCCW*weightVolumeCCW ) / ( weightVolumeCW + weightVolumeCCW )
      yVolWeighted <- ( yCartVolCW*weightVolumeCW + yCartVolCCW*weightVolumeCCW ) / ( weightVolumeCW + weightVolumeCCW )
      outVol[,,,nVols] <- atan2( yVolWeighted, xVolWeighted ) + pi
    }
  }
  write.AFNI( sprintf('%s_model_wAverage_10.nii.gz', participantName), brk=outVol, origin=CWmod$origin, orient=CWmod$orient, defhead=CWmod$NI_head )
  instr <- sprintf('compShiftPhase_model.sh %s_model_wAverage_10.nii.gz 180 %s_model_wAverage_10_correct.nii.gz 6', participantName, participantName )
  
  # clean up
  setwd( dirToCheck )
  print( getwd() )
  system('rm *zzz*')
  system( sprintf('rm %s', cwParticipants[nPart] ) )
  system( sprintf('rm %s', ccwParticipants[nPart] ) )
  system( sprintf('rm %s_CW_blur.nii.gz', participantName ) )
  system( sprintf('rm %s_CCW_blur.nii.gz', participantName ) )
  
}

