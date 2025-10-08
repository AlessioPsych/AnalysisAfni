rm( list=ls() );
mainDir <- '/analyse/Project0226/dataRepository'
targetDir <- '/analyse/Project0226/KastnerModel/data_KastnerClassic'
anatomiesDir_whole <- '/analyse/Project0226/KastnerModel/anatomies_KastnerClassic'
setwd( mainDir )

kastnerDir <- dir( pattern='_KastnerVest' )
#setwd( 'anatomies' )
anatomiesDir <- dir( anatomiesDir_whole )
#anatomiesDir <- anatomiesDir[c(1,2,3,5,7,8,10)]

data.frame( kastnerDir, anatomiesDir )

generalPurposeDir <- Sys.getenv( x='AFNI_TOOLBOXDIRGENERALPURPOSE' )
afniInstallDir <- Sys.getenv( x='AFNI_INSTALLDIR' ) 
source( sprintf('%s/scaleData.R', generalPurposeDir) )
source( sprintf('%s/AFNIio.R', afniInstallDir ) )
library( pracma )
library( abind )
library( circular )

for (nSubjects in 1:length(anatomiesDir)) { #subjId <- 5 length(anatomiesDir)
  
  subjId <- nSubjects
  
  ################
  # kastner data #
  ################
  #setwd( mainDir )
  #setwd( kastnerDir[ subjId ] ); getwd()
  #system('rm -R zzz*')
  
  # clean data on the anatomy folder
  setwd( anatomiesDir_whole )
  setwd( anatomiesDir[ subjId ] ); getwd()
  system('rm -R zzz*')
  
  # back to the original epi folder (needs to be created)
  setwd( mainDir )
  setwd( kastnerDir[ subjId ] ); getwd()
  
  #instr <- 'computeMeanTs.sh kastner_topUp/'; system( instr );
  #instr <- '3dTstat -prefix zzz_meanEpi4Coreg_kastner.nii.gz meanTs.nii'; system( instr )
  #system('rm meanTs.nii')
  
  #instr <- 'computeMeanTs.sh kastner_topUp_detrend/'; system( instr );
  #instr <- sprintf( 'mv meanTs.nii %s/%s/zzz_meanTs_kastner_topUp_dt.nii.gz', targetDir, kastnerDir[ subjId ] ); system( instr );
  instr <- sprintf('rm -R %s/%s', targetDir, kastnerDir[ subjId ]); system( instr )
  instr <- sprintf('mkdir %s/%s', targetDir, kastnerDir[ subjId ]); system( instr )
  instr <- sprintf( 'cp -R kastner_topUp_detrend/ %s/%s/kastner_topUp_detrend/', targetDir, kastnerDir[ subjId ] ); system( instr );
  instr <- sprintf( 'cp -R coregistration_kastner/ %s/%s/coregistration_kastner/', targetDir, kastnerDir[ subjId ] ); system( instr );
  instr <- sprintf( 'cp meanEpi4Coreg_kastner.nii.gz %s/%s/meanEpi4Coreg_kastner.nii.gz', targetDir, kastnerDir[ subjId ] ); system( instr );

  # move to target dir
  setwd( targetDir  )
  setwd( kastnerDir[ subjId ] ); getwd()
  
  #coregister
  setwd('coregistration_kastner')
  system('sh coregistrationScript')
  setwd( '..' ); getwd()
  
  # gray matter mask
  instr <- sprintf( '3dAllineate -source %s/%s/volumetricData_layering_depth.nii.gz -prefix zzz_greyVolumetric_kastner.nii.gz -1Dmatrix_apply coregistration_kastner/invMat.1D -master meanEpi4Coreg_kastner.nii.gz -final NN', anatomiesDir_whole, anatomiesDir[ subjId ] )
  system( instr )
  instr <- '3dcalc -a zzz_greyVolumetric_kastner.nii.gz -expx \u0027step(a)*within(a,0.01,0.99)\u0027 -prefix zzz_grey_kastner.nii.gz'; system( instr )
  
  #separate conditions
  setwd('kastner_topUp_detrend/')
  filesBrik <- dir(pattern='*.BRIK')
  filesHead <- dir(pattern='*.HEAD')
  
  baseInstrCW <- '3dMean -prefix CW_mean.nii.gz '
  for ( cwFiles in seq( 1, length(filesBrik), 2 ) ) {
    baseInstrCW <- sprintf('%s %s', baseInstrCW, filesBrik[cwFiles] )
  }
  print(baseInstrCW); system( baseInstrCW );
  instr <- 'mv CW_mean.nii.gz ../zzz_CW_mean.nii.gz'; system( instr );
  
  baseInstrCCW <- '3dMean -prefix CCW_mean.nii.gz '
  for ( ccwFiles in seq( 2, length(filesBrik), 2 ) ) {
    baseInstrCCW <- sprintf('%s %s', baseInstrCCW, filesBrik[ccwFiles] )
  }
  print(baseInstrCCW); system( baseInstrCCW );
  instr <- 'mv CCW_mean.nii.gz ../zzz_CCW_mean.nii.gz'; system( instr );
  
  setwd( '..' ); getwd()
  
  #instr <- '3dTcat -prefix zzz_CW_mean_delay.nii.gz zzz_CW_mean.nii.gz[2..161] zzz_CW_mean.nii.gz[0..1]'; system( instr )
  #instr <- '3dTcat -prefix zzz_CCW_mean_delay.nii.gz zzz_CCW_mean.nii.gz[2..161] zzz_CCW_mean.nii.gz[0..1]'; system( instr )
  #instr <- '3dcalc -a zzz_CCW_mean_delay.nii.gz -expr \u0027a*-1\u0027 -prefix zzz_CCW_mean_delay_inv.nii.gz'; system( instr )
  #instr <- '3dMean -prefix zzz_avg_kastner.nii.gz zzz_CW_mean_delay.nii.gz zzz_CCW_mean_delay_inv.nii.gz'; system( instr )

  #blur in mask
  #instr <- '3dBlurInMask -mask zzz_grey_kastner.nii.gz -FWHM 0.001 -prefix zzz_avg_kastner_blur.nii.gz -input zzz_avg_kastner.nii.gz'; system( instr )

  #smooth time
  #instr <- '3dTsmooth -prefix zzz_avg_kastner_blur_sm.nii.gz -lin zzz_avg_kastner_blur.nii.gz'; system( instr )

  #instr <- 'compTravellingWave.sh zzz_avg_kastner_blur_sm.nii.gz 5'; system( instr )
  #instr <- 'mv trWave.nii.gz zzz_trWave05.nii.gz'; system( instr )
  #instr <- 'mv periodogram.nii.gz zzz_periodogram05.nii.gz'; system( instr )

  #instr <- 'compTravellingWave.sh zzz_avg_kastner_blur_sm.nii.gz 10'; system( instr )
  #instr <- 'mv trWave.nii.gz zzz_trWave10.nii.gz'; system( instr )
  #instr <- 'mv periodogram.nii.gz zzz_periodogram10.nii.gz'; system( instr )
  
  #smooth time, CW, phase encoded analysis
  #instr <- '3dTsmooth -prefix zzz_CW_mean_sm.nii.gz -lin zzz_CW_mean.nii.gz'; system( instr )

  #blur in mask CW
  instr <- '3dBlurInMask -mask zzz_grey_kastner.nii.gz -FWHM 1.5 -prefix zzz_CW_mean_blur.nii.gz -input zzz_CW_mean.nii.gz'; system( instr )
  #blur in mask CCW
  instr <- '3dBlurInMask -mask zzz_grey_kastner.nii.gz -FWHM 1.5 -prefix zzz_CCW_mean_blur.nii.gz -input zzz_CCW_mean.nii.gz'; system( instr )
  
  #smooth time
  instr <- '3dTsmooth -prefix zzz_CW_mean_blur_sm.nii.gz -lin zzz_CW_mean_blur.nii.gz'; system( instr )
  instr <- '3dTsmooth -prefix zzz_CCW_mean_blur_sm.nii.gz -lin zzz_CCW_mean_blur.nii.gz'; system( instr )

  #travelling wave, 05 cycles, CW
  instr <- 'compTravellingWave.sh zzz_CW_mean_blur_sm.nii.gz 5'; system( instr )
  instr <- 'mv trWave.nii.gz zzz_trWave05_CW.nii.gz'; system( instr )
  instr <- 'mv periodogram.nii.gz zzz_periodogram05_CW.nii.gz'; system( instr )
  
  #travelling wave, 10 cycles, CW
  instr <- 'compTravellingWave.sh zzz_CW_mean_blur_sm.nii.gz 10'; system( instr )
  instr <- 'mv trWave.nii.gz zzz_trWave10_CW.nii.gz'; system( instr )
  instr <- 'mv periodogram.nii.gz zzz_periodogram10_CW.nii.gz'; system( instr )

  #travelling wave, 05 cycle, CCW
  instr <- 'compTravellingWave.sh zzz_CCW_mean_blur_sm.nii.gz 5'; system( instr )
  instr <- 'mv trWave.nii.gz zzz_trWave05_CCW.nii.gz'; system( instr )
  instr <- 'mv periodogram.nii.gz zzz_periodogram05_CCW.nii.gz'; system( instr )

  #travelling wave, 10 cycle, CCW
  instr <- 'compTravellingWave.sh zzz_CCW_mean_blur_sm.nii.gz 10'; system( instr )
  instr <- 'mv trWave.nii.gz zzz_trWave10_CCW.nii.gz'; system( instr )
  instr <- 'mv periodogram.nii.gz zzz_periodogram10_CCW.nii.gz'; system( instr )
  
  #shift phases from phase-spec coherence to average between CW and CCW (shift 90deg for haemodynamic delay, invert CW, 180deg to match with CCW)
  instr <- ' compShiftPhase_model.sh zzz_trWave05_CW.nii.gz 90 zzz_trWave05_CW_shift.nii.gz 3'; system( instr )
  instr <- ' compShiftPhase_model.sh zzz_trWave05_CCW.nii.gz 90 zzz_trWave05_CCW_shift.nii.gz 3'; system( instr )
  instr <- ' compShiftPhase_model.sh zzz_trWave05_CW_shift.nii.gz 180 zzz_trWave05_CW_shift_register.nii.gz 3'; system( instr )

  instr <- ' compShiftPhase_model.sh zzz_trWave10_CW.nii.gz 90 zzz_trWave10_CW_shift.nii.gz 3'; system( instr )
  instr <- ' compShiftPhase_model.sh zzz_trWave10_CCW.nii.gz 90 zzz_trWave10_CCW_shift.nii.gz 3'; system( instr )
  instr <- ' compShiftPhase_model.sh zzz_trWave10_CW_shift.nii.gz 180 zzz_trWave10_CW_shift_register.nii.gz 3'; system( instr )
  
  # (weighted) average phase-encoded, shifted results, 05 cycles:
  CWmod <- read.AFNI('zzz_trWave05_CW_shift_register.nii.gz')
  CCWmod <- read.AFNI('zzz_trWave05_CCW_shift.nii.gz')
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
      xCartVolCW <- cos( CWmodVol[,,,nVols] ) # convert to cartesian assuming a radius or 1
      yCartVolCW <- sin( CWmodVol[,,,nVols] ) # convert to cartesian assuming a radius or 1
      xCartVolCCW <- cos( CCWmodVol[,,,nVols] ) # convert to cartesian assuming a radius or 1
      yCartVolCCW <- sin( CCWmodVol[,,,nVols] ) # convert to cartesian assuming a radius or 1
      xVolWeighted <- ( xCartVolCW*weightVolumeCW + xCartVolCCW*weightVolumeCCW ) / ( weightVolumeCW + weightVolumeCCW )
      yVolWeighted <- ( yCartVolCW*weightVolumeCW + yCartVolCCW*weightVolumeCCW ) / ( weightVolumeCW + weightVolumeCCW )
      outVol[,,,nVols] <- atan2( yVolWeighted, xVolWeighted ) + pi
    }
  }
  write.AFNI( 'zzz_phase_encoded_wAverage_05.nii.gz', brk=outVol, origin=CWmod$origin, orient=CWmod$orient, defhead=CWmod$NI_head )
  # shift output, atan2 gives the inverted output
  instr <- ' compShiftPhase_model.sh zzz_phase_encoded_wAverage_05.nii.gz 180 zzz_phase_encoded_wAverage_05_correct.nii.gz 3'; system( instr )

  # (weighted) average phase-encoded, shifted results, 10 cycles:
  CWmod <- read.AFNI('zzz_trWave10_CW_shift_register.nii.gz')
  CCWmod <- read.AFNI('zzz_trWave10_CCW_shift.nii.gz')
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
      xCartVolCW <- cos( CWmodVol[,,,nVols] ) # convert to cartesian assuming a radius or 1
      yCartVolCW <- sin( CWmodVol[,,,nVols] ) # convert to cartesian assuming a radius or 1
      xCartVolCCW <- cos( CCWmodVol[,,,nVols] ) # convert to cartesian assuming a radius or 1
      yCartVolCCW <- sin( CCWmodVol[,,,nVols] ) # convert to cartesian assuming a radius or 1
      xVolWeighted <- ( xCartVolCW*weightVolumeCW + xCartVolCCW*weightVolumeCCW ) / ( weightVolumeCW + weightVolumeCCW )
      yVolWeighted <- ( yCartVolCW*weightVolumeCW + yCartVolCCW*weightVolumeCCW ) / ( weightVolumeCW + weightVolumeCCW )
      outVol[,,,nVols] <- atan2( yVolWeighted, xVolWeighted ) + pi
    }
  }
  write.AFNI( 'zzz_phase_encoded_wAverage_10.nii.gz', brk=outVol, origin=CWmod$origin, orient=CWmod$orient, defhead=CWmod$NI_head )
  # shift output, atan2 gives the inverted output
  instr <- ' compShiftPhase_model.sh zzz_phase_encoded_wAverage_10.nii.gz 180 zzz_phase_encoded_wAverage_10_correct.nii.gz 3'; system( instr )
  
  # #kastner model prediction, sinewave
  # instr <- 'Rscript /analyse/Project0226/dataToSort/eyeMovPoster/codeFit/generatePredictionsKastnerSinewave_PSOCK.R zzz_CW_mean_blur_sm.nii.gz zzz_kastner_sinewave/ 0 cw 5'; system( instr )
  #kastner model prediction CW
  #instr <- 'Rscript /analyse/Project0226/dataToSort/eyeMovPoster/codeFit/generatePredictionsKastner_PSOCK.R zzz_CW_mean_blur_sm.nii.gz zzz_kastnerCW/ 0 cw'; system( instr )
  #kastner model prediction CCW
  #instr <- 'Rscript /analyse/Project0226/dataToSort/eyeMovPoster/codeFit/generatePredictionsKastner_PSOCK.R zzz_CCW_mean_blur_sm.nii.gz zzz_kastnerCCW/ 0 ccw'; system( instr )
  #kastner model prediction whole
  #instr <- 'Rscript /analyse/Project0226/dataToSort/eyeMovPoster/codeFit/generatePredictionsKastner_whole.R zzz_whole_blur.nii.gz zzz_kastner_whole/ 0'; system( instr )
  
  #kastner model prediction 2 components CW
  instr <- 'Rscript ../../code/generatePredictionsKastner_2Components_differentDirection_PSOCK.R zzz_CW_mean_blur_sm.nii.gz zzz_kastnerCW_2components_both/ 0 cw both'; system( instr )
  #kastner model prediction 2 components CCW
  instr <- 'Rscript ../../code/generatePredictionsKastner_2Components_differentDirection_PSOCK.R zzz_CCW_mean_blur_sm.nii.gz zzz_kastnerCCW_2components_both/ 0 ccw both'; system( instr )
  
  #kastner model prediction 2 components CW 10
  instr <- 'Rscript ../../code/generatePredictionsKastner_2Components_differentDirection_PSOCK.R zzz_CW_mean_blur_sm.nii.gz zzz_kastnerCW_2components_10/ 0 cw 10'; system( instr )
  #kastner model prediction 2 components CCW 10
  instr <- 'Rscript ../../code/generatePredictionsKastner_2Components_differentDirection_PSOCK.R zzz_CCW_mean_blur_sm.nii.gz zzz_kastnerCCW_2components_10/ 0 ccw 10'; system( instr )

  #kastner model prediction 2 components CW 5
  instr <- 'Rscript ../../code/generatePredictionsKastner_2Components_differentDirection_PSOCK.R zzz_CW_mean_blur_sm.nii.gz zzz_kastnerCW_2components_05/ 0 cw 5'; system( instr )
  #kastner model prediction 2 components CCW 5
  instr <- 'Rscript ../../code/generatePredictionsKastner_2Components_differentDirection_PSOCK.R zzz_CCW_mean_blur_sm.nii.gz zzz_kastnerCCW_2components_05/ 0 ccw 5'; system( instr )
  
  # mask blur
  instr <- '3dTstat -prefix zzz_mean_blur.nii.gz zzz_CW_mean_blur.nii.gz'; system( instr )
  instr <- sprintf('3dcalc -a zzz_mean_blur.nii.gz -expr \u0027step(abs(a))\u0027 -prefix zzz_mask_blur.nii.gz'); system( instr )
  
  #kastner model, CW
  #instr <- 'Rscript ../../code/fitEyeKinematics_PSOCK.R zzz_CW_mean_blur_sm.nii.gz zzz_mask_blur.nii.gz zzz_kastnerCW/ zzz_CW'; system( instr )
  #kastner model, CCW
  #instr <- 'Rscript ../../code/fitEyeKinematics_PSOCK.R zzz_CCW_mean_blur_sm.nii.gz zzz_mask_blur.nii.gz zzz_kastnerCCW/ zzz_CCW'; system( instr )
  # #kastner model, CW, sinewave
  # instr <- 'Rscript /analyse/Project0226/dataToSort/eyeMovPoster/codeFit/fitEyeKinematics_PSOCK.R zzz_CW_mean_blur_sm.nii.gz zzz_mask_blur.nii.gz zzz_kastner_sinewave/ zzz_CW_sinewave'; system( instr )
  # #kastner model, CCW, sinewave
  # instr <- 'Rscript /analyse/Project0226/dataToSort/eyeMovPoster/codeFit/fitEyeKinematics_PSOCK.R zzz_CCW_mean_blur_sm.nii.gz zzz_mask_blur.nii.gz zzz_kastner_sinewave/ zzz_CCW_sinewave'; system( instr )
  
  #kastner model 2 components, CW both
  instr <- 'Rscript ../../code/fitEyeKinematics_PSOCK.R zzz_CW_mean_blur_sm.nii.gz zzz_mask_blur.nii.gz zzz_kastnerCW_2components_both/ zzz_CW_2components_both kinematicModel'; system( instr )
  #kastner model 2 components, CCW both
  instr <- 'Rscript ../../code/fitEyeKinematics_PSOCK.R zzz_CCW_mean_blur_sm.nii.gz zzz_mask_blur.nii.gz zzz_kastnerCCW_2components_both/ zzz_CCW_2components_both kinematicModel'; system( instr )

  #kastner model 2 components, CW 10
  instr <- 'Rscript ../../code/fitEyeKinematics_PSOCK.R zzz_CW_mean_blur_sm.nii.gz zzz_mask_blur.nii.gz zzz_kastnerCW_2components_10/ zzz_CW_2components_10 kinematicModel'; system( instr )
  #kastner model 2 components, CCW 10
  instr <- 'Rscript ../../code/fitEyeKinematics_PSOCK.R zzz_CCW_mean_blur_sm.nii.gz zzz_mask_blur.nii.gz zzz_kastnerCCW_2components_10/ zzz_CCW_2components_10 kinematicModel'; system( instr )

  #kastner model 2 components, CW 5
  instr <- 'Rscript ../../code/fitEyeKinematics_PSOCK.R zzz_CW_mean_blur_sm.nii.gz zzz_mask_blur.nii.gz zzz_kastnerCW_2components_05/ zzz_CW_2components_05 kinematicModel'; system( instr )
  #kastner model 2 components, CCW 5
  instr <- 'Rscript ../../code/fitEyeKinematics_PSOCK.R zzz_CCW_mean_blur_sm.nii.gz zzz_mask_blur.nii.gz zzz_kastnerCCW_2components_05/ zzz_CCW_2components_05 kinematicModel'; system( instr )

  
  # (weighted) average modelling results, two components:
  CWmod <- read.AFNI('zzz_CW_2components_both_params.nii.gz')
  CCWmod <- read.AFNI('zzz_CCW_2components_both_params.nii.gz')
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
      xCartVolCW <- cos( CWmodVol[,,,nVols] ) # convert to cartesian assuming a radius or 1
      yCartVolCW <- sin( CWmodVol[,,,nVols] ) # convert to cartesian assuming a radius or 1
      xCartVolCCW <- cos( CCWmodVol[,,,nVols] ) # convert to cartesian assuming a radius or 1
      yCartVolCCW <- sin( CCWmodVol[,,,nVols] ) # convert to cartesian assuming a radius or 1
      xVolWeighted <- ( xCartVolCW*weightVolumeCW + xCartVolCCW*weightVolumeCCW ) / ( weightVolumeCW + weightVolumeCCW )
      yVolWeighted <- ( yCartVolCW*weightVolumeCW + yCartVolCCW*weightVolumeCCW ) / ( weightVolumeCW + weightVolumeCCW )
      outVol[,,,nVols] <- atan2( yVolWeighted, xVolWeighted ) + pi
    }
  }
  write.AFNI( 'zzz_model_wAverage_both.nii.gz', brk=outVol, origin=CWmod$origin, orient=CWmod$orient, defhead=CWmod$NI_head )
  instr <- ' compShiftPhase_model.sh zzz_model_wAverage_both.nii.gz 180 zzz_model_wAverage_both_correct.nii.gz 6'; system( instr )
  
  # (weighted) average modelling results, 05 cycles:
  CWmod <- read.AFNI('zzz_CW_2components_05_params.nii.gz')
  CCWmod <- read.AFNI('zzz_CCW_2components_05_params.nii.gz')
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
      xCartVolCW <- cos( CWmodVol[,,,nVols] ) # convert to cartesian assuming a radius or 1
      yCartVolCW <- sin( CWmodVol[,,,nVols] ) # convert to cartesian assuming a radius or 1
      xCartVolCCW <- cos( CCWmodVol[,,,nVols] ) # convert to cartesian assuming a radius or 1
      yCartVolCCW <- sin( CCWmodVol[,,,nVols] ) # convert to cartesian assuming a radius or 1
      xVolWeighted <- ( xCartVolCW*weightVolumeCW + xCartVolCCW*weightVolumeCCW ) / ( weightVolumeCW + weightVolumeCCW )
      yVolWeighted <- ( yCartVolCW*weightVolumeCW + yCartVolCCW*weightVolumeCCW ) / ( weightVolumeCW + weightVolumeCCW )
      outVol[,,,nVols] <- atan2( yVolWeighted, xVolWeighted ) + pi
    }
  }
  write.AFNI( 'zzz_model_wAverage_05.nii.gz', brk=outVol, origin=CWmod$origin, orient=CWmod$orient, defhead=CWmod$NI_head )
  instr <- ' compShiftPhase_model.sh zzz_model_wAverage_05.nii.gz 180 zzz_model_wAverage_05_correct.nii.gz 6'; system( instr )

  # (weighted) average modelling results, 10 cycles:
  CWmod <- read.AFNI('zzz_CW_2components_10_params.nii.gz')
  CCWmod <- read.AFNI('zzz_CCW_2components_10_params.nii.gz')
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
      xCartVolCW <- cos( CWmodVol[,,,nVols] ) # convert to cartesian assuming a radius or 1
      yCartVolCW <- sin( CWmodVol[,,,nVols] ) # convert to cartesian assuming a radius or 1
      xCartVolCCW <- cos( CCWmodVol[,,,nVols] ) # convert to cartesian assuming a radius or 1
      yCartVolCCW <- sin( CCWmodVol[,,,nVols] ) # convert to cartesian assuming a radius or 1
      xVolWeighted <- ( xCartVolCW*weightVolumeCW + xCartVolCCW*weightVolumeCCW ) / ( weightVolumeCW + weightVolumeCCW )
      yVolWeighted <- ( yCartVolCW*weightVolumeCW + yCartVolCCW*weightVolumeCCW ) / ( weightVolumeCW + weightVolumeCCW )
      outVol[,,,nVols] <- atan2( yVolWeighted, xVolWeighted ) + pi
    }
  }
  write.AFNI( 'zzz_model_wAverage_10.nii.gz', brk=outVol, origin=CWmod$origin, orient=CWmod$orient, defhead=CWmod$NI_head )
  instr <- ' compShiftPhase_model.sh zzz_model_wAverage_10.nii.gz 180 zzz_model_wAverage_10_correct.nii.gz 6'; system( instr )
  
  # copy results in the anatomy, modelling
  instr <- sprintf( '3dAllineate -source zzz_CCW_2components_both_params.nii.gz -master ../../anatomies_KastnerClassic/%s/anatCopy.nii.gz -prefix ../../anatomies_KastnerClassic/%s/zzz_kastner_model_both_params_anat_ccw.nii.gz -1Dmatrix_apply coregistration_kastner/coregMat.1D -final NN', anatomiesDir[ subjId ], anatomiesDir[ subjId ] )
  system( instr )
  instr <- sprintf( '3dAllineate -source zzz_CCW_2components_10_params.nii.gz -master ../../anatomies_KastnerClassic/%s/anatCopy.nii.gz -prefix ../../anatomies_KastnerClassic/%s/zzz_kastner_model_10_params_anat_ccw.nii.gz -1Dmatrix_apply coregistration_kastner/coregMat.1D -final NN', anatomiesDir[ subjId ], anatomiesDir[ subjId ] )
  system( instr )
  instr <- sprintf( '3dAllineate -source zzz_CCW_2components_05_params.nii.gz -master ../../anatomies_KastnerClassic/%s/anatCopy.nii.gz -prefix ../../anatomies_KastnerClassic/%s/zzz_kastner_model_05_params_anat_ccw.nii.gz -1Dmatrix_apply coregistration_kastner/coregMat.1D -final NN', anatomiesDir[ subjId ], anatomiesDir[ subjId ] )
  system( instr )
  instr <- sprintf( '3dAllineate -source zzz_CW_2components_both_params.nii.gz -master ../../anatomies_KastnerClassic/%s/anatCopy.nii.gz -prefix ../../anatomies_KastnerClassic/%s/zzz_kastner_model_both_params_anat_cw.nii.gz -1Dmatrix_apply coregistration_kastner/coregMat.1D -final NN', anatomiesDir[ subjId ], anatomiesDir[ subjId ] )
  system( instr )
  instr <- sprintf( '3dAllineate -source zzz_CW_2components_10_params.nii.gz -master ../../anatomies_KastnerClassic/%s/anatCopy.nii.gz -prefix ../../anatomies_KastnerClassic/%s/zzz_kastner_model_10_params_anat_cw.nii.gz -1Dmatrix_apply coregistration_kastner/coregMat.1D -final NN', anatomiesDir[ subjId ], anatomiesDir[ subjId ] )
  system( instr )
  instr <- sprintf( '3dAllineate -source zzz_CW_2components_05_params.nii.gz -master ../../anatomies_KastnerClassic/%s/anatCopy.nii.gz -prefix ../../anatomies_KastnerClassic/%s/zzz_kastner_model_05_params_anat_cw.nii.gz -1Dmatrix_apply coregistration_kastner/coregMat.1D -final NN', anatomiesDir[ subjId ], anatomiesDir[ subjId ] )
  system( instr )
  instr <- sprintf( '3dAllineate -source zzz_model_wAverage_05_correct.nii.gz -master ../../anatomies_KastnerClassic/%s/anatCopy.nii.gz -prefix ../../anatomies_KastnerClassic/%s/zzz_model_wAverage_05_correct_anat.nii.gz -1Dmatrix_apply coregistration_kastner/coregMat.1D -final NN', anatomiesDir[ subjId ], anatomiesDir[ subjId ] )
  system( instr )
  instr <- sprintf( '3dAllineate -source zzz_model_wAverage_10_correct.nii.gz -master ../../anatomies_KastnerClassic/%s/anatCopy.nii.gz -prefix ../../anatomies_KastnerClassic/%s/zzz_model_wAverage_10_correct_anat.nii.gz -1Dmatrix_apply coregistration_kastner/coregMat.1D -final NN', anatomiesDir[ subjId ], anatomiesDir[ subjId ] )
  system( instr )
  instr <- sprintf( '3dAllineate -source zzz_model_wAverage_both_correct.nii.gz -master ../../anatomies_KastnerClassic/%s/anatCopy.nii.gz -prefix ../../anatomies_KastnerClassic/%s/zzz_model_wAverage_both_correct_anat.nii.gz -1Dmatrix_apply coregistration_kastner/coregMat.1D -final NN', anatomiesDir[ subjId ], anatomiesDir[ subjId ] )
  system( instr )
  
  # copy results in the anatomy, phase-encoded
  instr <- sprintf( '3dAllineate -source zzz_trWave05_CW.nii.gz -master ../../anatomies_KastnerClassic/%s/anatCopy.nii.gz -prefix ../../anatomies_KastnerClassic/%s/zzz_trWave05_CW_anat.nii.gz -1Dmatrix_apply coregistration_kastner/coregMat.1D -final NN', anatomiesDir[ subjId ], anatomiesDir[ subjId ] )
  system( instr )
  instr <- sprintf( '3dAllineate -source zzz_trWave05_CCW.nii.gz -master ../../anatomies_KastnerClassic/%s/anatCopy.nii.gz -prefix ../../anatomies_KastnerClassic/%s/zzz_trWave05_CCW_anat.nii.gz -1Dmatrix_apply coregistration_kastner/coregMat.1D -final NN', anatomiesDir[ subjId ], anatomiesDir[ subjId ] )
  system( instr )
  instr <- sprintf( '3dAllineate -source zzz_phase_encoded_wAverage_05_correct.nii.gz -master ../../anatomies_KastnerClassic/%s/anatCopy.nii.gz -prefix ../../anatomies_KastnerClassic/%s/zzz_phase_encoded_wAverage_05_correct_anat.nii.gz -1Dmatrix_apply coregistration_kastner/coregMat.1D -final NN', anatomiesDir[ subjId ], anatomiesDir[ subjId ] )
  system( instr )
  instr <- sprintf( '3dAllineate -source zzz_trWave10_CW.nii.gz -master ../../anatomies_KastnerClassic/%s/anatCopy.nii.gz -prefix ../../anatomies_KastnerClassic/%s/zzz_trWave10_CW_anat.nii.gz -1Dmatrix_apply coregistration_kastner/coregMat.1D -final NN', anatomiesDir[ subjId ], anatomiesDir[ subjId ] )
  system( instr )
  instr <- sprintf( '3dAllineate -source zzz_trWave10_CCW.nii.gz -master ../../anatomies_KastnerClassic/%s/anatCopy.nii.gz -prefix ../../anatomies_KastnerClassic/%s/zzz_trWave10_CCW_anat.nii.gz -1Dmatrix_apply coregistration_kastner/coregMat.1D -final NN', anatomiesDir[ subjId ], anatomiesDir[ subjId ] )
  system( instr )
  instr <- sprintf( '3dAllineate -source zzz_phase_encoded_wAverage_10_correct.nii.gz -master ../../anatomies_KastnerClassic/%s/anatCopy.nii.gz -prefix ../../anatomies_KastnerClassic/%s/zzz_phase_encoded_wAverage_10_correct_anat.nii.gz -1Dmatrix_apply coregistration_kastner/coregMat.1D -final NN', anatomiesDir[ subjId ], anatomiesDir[ subjId ] )
  system( instr )
  
  # move to anatomical directory and interpolate on surface
  setwd( sprintf('../../anatomies_KastnerClassic/%s/', anatomiesDir[ subjId ]) ); getwd()
  instr <- sprintf( '3dcalc -a zzz_kastner_model_both_params_anat_ccw.nii.gz -exp \u0027a+0.0001\u0027 -prefix zzz_kastner_model_both_params_anat_ccw_add.nii.gz'); system( instr )  
  instr <- sprintf( '3dcalc -a zzz_kastner_model_10_params_anat_ccw.nii.gz -exp \u0027a+0.0001\u0027 -prefix zzz_kastner_model_10_params_anat_ccw_add.nii.gz'); system( instr )  
  instr <- sprintf( '3dcalc -a zzz_kastner_model_05_params_anat_ccw.nii.gz -exp \u0027a+0.0001\u0027 -prefix zzz_kastner_model_05_params_anat_ccw_add.nii.gz'); system( instr )  
  
  instr <- sprintf( '3dcalc -a zzz_kastner_model_both_params_anat_cw.nii.gz -exp \u0027a+0.0001\u0027 -prefix zzz_kastner_model_both_params_anat_cw_add.nii.gz'); system( instr )  
  instr <- sprintf( '3dcalc -a zzz_kastner_model_10_params_anat_cw.nii.gz -exp \u0027a+0.0001\u0027 -prefix zzz_kastner_model_10_params_anat_cw_add.nii.gz'); system( instr )  
  instr <- sprintf( '3dcalc -a zzz_kastner_model_05_params_anat_cw.nii.gz -exp \u0027a+0.0001\u0027 -prefix zzz_kastner_model_05_params_anat_cw_add.nii.gz'); system( instr )  

  instr <- sprintf( '3dcalc -a zzz_model_wAverage_05_correct_anat.nii.gz -exp \u0027a+0.0001\u0027 -prefix zzz_model_wAverage_05_correct_anat_add.nii.gz'); system( instr )  
  instr <- sprintf( '3dcalc -a zzz_model_wAverage_10_correct_anat.nii.gz -exp \u0027a+0.0001\u0027 -prefix zzz_model_wAverage_10_correct_anat_add.nii.gz'); system( instr )  
  instr <- sprintf( '3dcalc -a zzz_model_wAverage_both_correct_anat.nii.gz -exp \u0027a+0.0001\u0027 -prefix zzz_model_wAverage_both_correct_anat_add.nii.gz'); system( instr )  
  
  instr <- sprintf( '3dcalc -a zzz_trWave05_CW_anat.nii.gz -exp \u0027a+0.0001\u0027 -prefix zzz_trWave05_CW_anat_add.nii.gz'); system( instr )  
  instr <- sprintf( '3dcalc -a zzz_trWave05_CCW_anat.nii.gz -exp \u0027a+0.0001\u0027 -prefix zzz_trWave05_CCW_anat_add.nii.gz'); system( instr )  
  instr <- sprintf( '3dcalc -a zzz_trWave10_CW_anat.nii.gz -exp \u0027a+0.0001\u0027 -prefix zzz_trWave10_CW_anat_add.nii.gz'); system( instr )  
  instr <- sprintf( '3dcalc -a zzz_trWave10_CCW_anat.nii.gz -exp \u0027a+0.0001\u0027 -prefix zzz_trWave10_CCW_anat_add.nii.gz'); system( instr )  

  instr <- sprintf( '3dcalc -a zzz_phase_encoded_wAverage_05_correct_anat.nii.gz -exp \u0027a+0.0001\u0027 -prefix zzz_phase_encoded_wAverage_05_correct_anat_add.nii.gz'); system( instr )  
  instr <- sprintf( '3dcalc -a zzz_phase_encoded_wAverage_10_correct_anat.nii.gz -exp \u0027a+0.0001\u0027 -prefix zzz_phase_encoded_wAverage_10_correct_anat_add.nii.gz'); system( instr )  
  instr <- sprintf( '3dcalc -a zzz_phase_encoded_wAverage_both_correct_anat.nii.gz -exp \u0027a+0.0001\u0027 -prefix zzz_phase_encoded_wAverage_both_correct_anat_add.nii.gz'); system( instr )  
  
  #instr <- sprintf( 'interpolateVolumeOnSurface.sh zzz_kastner_model_both_params_anat_ccw_add.nii.gz surfaces_folder_left/'  ); system( instr )
  #instr <- sprintf( 'interpolateVolumeOnSurface.sh zzz_kastner_model_both_params_anat_ccw_add.nii.gz surfaces_folder_right/'  ); system( instr )
  instr <- sprintf( 'interpolateVolumeOnSurface.sh zzz_kastner_model_10_params_anat_ccw_add.nii.gz surfaces_folder_left/'  ); system( instr )
  instr <- sprintf( 'interpolateVolumeOnSurface.sh zzz_kastner_model_10_params_anat_ccw_add.nii.gz surfaces_folder_right/'  ); system( instr )
  instr <- sprintf( 'interpolateVolumeOnSurface.sh zzz_kastner_model_05_params_anat_ccw_add.nii.gz surfaces_folder_left/'  ); system( instr )
  instr <- sprintf( 'interpolateVolumeOnSurface.sh zzz_kastner_model_05_params_anat_ccw_add.nii.gz surfaces_folder_right/'  ); system( instr )

  #instr <- sprintf( 'interpolateVolumeOnSurface.sh zzz_kastner_model_both_params_anat_cw_add.nii.gz surfaces_folder_left/'  ); system( instr )
  #instr <- sprintf( 'interpolateVolumeOnSurface.sh zzz_kastner_model_both_params_anat_cw_add.nii.gz surfaces_folder_right/'  ); system( instr )
  instr <- sprintf( 'interpolateVolumeOnSurface.sh zzz_kastner_model_10_params_anat_cw_add.nii.gz surfaces_folder_left/'  ); system( instr )
  instr <- sprintf( 'interpolateVolumeOnSurface.sh zzz_kastner_model_10_params_anat_cw_add.nii.gz surfaces_folder_right/'  ); system( instr )
  instr <- sprintf( 'interpolateVolumeOnSurface.sh zzz_kastner_model_05_params_anat_cw_add.nii.gz surfaces_folder_left/'  ); system( instr )
  instr <- sprintf( 'interpolateVolumeOnSurface.sh zzz_kastner_model_05_params_anat_cw_add.nii.gz surfaces_folder_right/'  ); system( instr )

  instr <- sprintf( 'interpolateVolumeOnSurface.sh zzz_model_wAverage_05_correct_anat_add.nii.gz surfaces_folder_left/'  ); system( instr )
  instr <- sprintf( 'interpolateVolumeOnSurface.sh zzz_model_wAverage_05_correct_anat_add.nii.gz surfaces_folder_right/'  ); system( instr )
  instr <- sprintf( 'interpolateVolumeOnSurface.sh zzz_model_wAverage_10_correct_anat_add.nii.gz surfaces_folder_left/'  ); system( instr )
  instr <- sprintf( 'interpolateVolumeOnSurface.sh zzz_model_wAverage_10_correct_anat_add.nii.gz surfaces_folder_right/'  ); system( instr )
  instr <- sprintf( 'interpolateVolumeOnSurface.sh zzz_model_wAverage_both_correct_anat_add.nii.gz surfaces_folder_left/'  ); system( instr )
  instr <- sprintf( 'interpolateVolumeOnSurface.sh zzz_model_wAverage_both_correct_anat_add.nii.gz surfaces_folder_right/'  ); system( instr )
  
  instr <- sprintf( 'interpolateVolumeOnSurface.sh zzz_trWave05_CW_anat_add.nii.gz surfaces_folder_left/'  ); system( instr )
  instr <- sprintf( 'interpolateVolumeOnSurface.sh zzz_trWave05_CW_anat_add.nii.gz surfaces_folder_right/'  ); system( instr )
  instr <- sprintf( 'interpolateVolumeOnSurface.sh zzz_trWave05_CCW_anat_add.nii.gz surfaces_folder_left/'  ); system( instr )
  instr <- sprintf( 'interpolateVolumeOnSurface.sh zzz_trWave05_CCW_anat_add.nii.gz surfaces_folder_right/'  ); system( instr )
  instr <- sprintf( 'interpolateVolumeOnSurface.sh zzz_trWave10_CW_anat_add.nii.gz surfaces_folder_left/'  ); system( instr )
  instr <- sprintf( 'interpolateVolumeOnSurface.sh zzz_trWave10_CW_anat_add.nii.gz surfaces_folder_right/'  ); system( instr )
  instr <- sprintf( 'interpolateVolumeOnSurface.sh zzz_trWave10_CCW_anat_add.nii.gz surfaces_folder_left/'  ); system( instr )
  instr <- sprintf( 'interpolateVolumeOnSurface.sh zzz_trWave10_CCW_anat_add.nii.gz surfaces_folder_right/'  ); system( instr )

  instr <- sprintf( 'interpolateVolumeOnSurface.sh zzz_phase_encoded_wAverage_05_correct_anat_add.nii.gz surfaces_folder_left/'  ); system( instr )
  instr <- sprintf( 'interpolateVolumeOnSurface.sh zzz_phase_encoded_wAverage_05_correct_anat_add.nii.gz surfaces_folder_right/'  ); system( instr )
  instr <- sprintf( 'interpolateVolumeOnSurface.sh zzz_phase_encoded_wAverage_10_correct_anat_add.nii.gz surfaces_folder_left/'  ); system( instr )
  instr <- sprintf( 'interpolateVolumeOnSurface.sh zzz_phase_encoded_wAverage_10_correct_anat_add.nii.gz surfaces_folder_right/'  ); system( instr )
  
  # instr <- 'compShiftPhase_model.sh zzz_CW_sinewave_params.nii.gz 33 zzz_CW_sinewave_params_shft33.nii.gz 6'; system( instr )
  # instr <- 'compShiftPhase_model.sh zzz_CCW_sinewave_params.nii.gz 0 zzz_CCW_sinewave_params_shft33.nii.gz 6'; system( instr )
  # instr <- 'compShiftPhase_model.sh zzz_CCW_sinewave_params_shft33.nii.gz 33 zzz_CCW_sinewave_params_shft33_inv.nii.gz 6'; system( instr )
  
  #set.seed(123); 
  #par(mfrow=c(1,3)); 
  #phv <- rnorm(10,0.5,0.9); 
  #hist( phv, breaks=9, xlim=c(-3.14, 3.14) ); 
  #phvt01 <- (phv+pi+0) ; 
  #hist( phvt01, breaks=9, xlim=c(0, 6.28) );
  #phvt02 <- (phv+0) %% (2*pi) ; 
  #hist( phvt02, breaks=9, xlim=c(0, 6.28) );
  
  #kastner model, whole
  #instr <- 'Rscript /analyse/Project0226/dataToSort/eyeMovPoster/codeFit/fitEyeKinematics.R zzz_whole_blur_smooth.nii.gz zzz_mask_blur.nii.gz zzz_kastner_whole/ zzz_whole_smooth'; system( instr )
  
  # CWmod <- read.AFNI('zzz_CW_params.nii.gz')
  # CCWmod <- read.AFNI('zzz_CCW_params.nii.gz')
  # CWmodVol <- CWmod$brk
  # CCWmodVol <- CCWmod$brk
  # outVol <- array(0, dim( CCWmod$brk ) )
  # for (nVols in 1:dim( CCWmod$brk )[4] ) {
  #   if (nVols!=6) { outVol[,,,nVols] <- ( CWmodVol[,,,nVols] + CCWmodVol[,,,nVols] ) / 2 }
  #   if (nVols==6) { outVol[,,,nVols] <- atan2( ( sin( CWmodVol[,,,nVols] ) + sin( CCWmodVol[,,,nVols] ) ) / 2, 
  #                                              ( cos( CWmodVol[,,,nVols] ) + cos( CCWmodVol[,,,nVols] ) ) / 2 )    }
  # }
  # write.AFNI( 'zzz_model_average.nii.gz', brk=outVol, origin=CWmod$origin, orient=CWmod$orient, defhead=CWmod$NI_head )
  # 
  # CWmod <- read.AFNI('zzz_CW_2components_05_params.nii.gz')
  # CCWmod <- read.AFNI('zzz_CCW_2components_05_params.nii.gz')
  # CWmodVol <- CWmod$brk
  # CCWmodVol <- CCWmod$brk
  # outVol <- array(0, dim( CCWmod$brk ) )
  # for (nVols in 1:dim( CCWmod$brk )[4] ) {
  #  if (nVols!=6) { outVol[,,,nVols] <- ( CWmodVol[,,,nVols] + CCWmodVol[,,,nVols] ) / 2 }
  #  if (nVols==6) { outVol[,,,nVols] <- atan2( ( sin( CWmodVol[,,,nVols] ) + sin( CCWmodVol[,,,nVols] ) ) / 2, 
  #                                             ( cos( CWmodVol[,,,nVols] ) + cos( CCWmodVol[,,,nVols] ) ) / 2 )    }
  # }
  # write.AFNI( 'zzz_model_average_2components_05.nii.gz', brk=outVol, origin=CWmod$origin, orient=CWmod$orient, defhead=CWmod$NI_head )
  # 
  # CWmod <- read.AFNI('zzz_CW_2components_10_params.nii.gz')
  # CCWmod <- read.AFNI('zzz_CCW_2components_10_params.nii.gz')
  # CWmodVol <- CWmod$brk
  # CCWmodVol <- CCWmod$brk
  # outVol <- array(0, dim( CCWmod$brk ) )
  # for (nVols in 1:dim( CCWmod$brk )[4] ) {
  #   if (nVols!=6) { outVol[,,,nVols] <- ( CWmodVol[,,,nVols] + CCWmodVol[,,,nVols] ) / 2 }
  #   if (nVols==6) { outVol[,,,nVols] <- atan2( ( sin( CWmodVol[,,,nVols] ) + sin( CCWmodVol[,,,nVols] ) ) / 2, 
  #                                              ( cos( CWmodVol[,,,nVols] ) + cos( CCWmodVol[,,,nVols] ) ) / 2 )    }
  # }
  # write.AFNI( 'zzz_model_average_2components_10.nii.gz', brk=outVol, origin=CWmod$origin, orient=CWmod$orient, defhead=CWmod$NI_head )
  # 
  # CWmod <- read.AFNI('zzz_CW_2components_both_params.nii.gz')
  # CCWmod <- read.AFNI('zzz_CCW_2components_both_params.nii.gz')
  # CWmodVol <- CWmod$brk
  # CCWmodVol <- CCWmod$brk
  # outVol <- array(0, dim( CCWmod$brk ) )
  # for (nVols in 1:dim( CCWmod$brk )[4] ) {
  #   if (nVols!=6) { outVol[,,,nVols] <- ( CWmodVol[,,,nVols] + CCWmodVol[,,,nVols] ) / 2 }
  #   if (nVols==6) { outVol[,,,nVols] <- atan2( ( sin( CWmodVol[,,,nVols] ) + sin( CCWmodVol[,,,nVols] ) ) / 2, 
  #                                              ( cos( CWmodVol[,,,nVols] ) + cos( CCWmodVol[,,,nVols] ) ) / 2 )    }
  # }
  # write.AFNI( 'zzz_model_average_2components_both.nii.gz', brk=outVol, origin=CWmod$origin, orient=CWmod$orient, defhead=CWmod$NI_head )
  # 
  

}

# simple model with one covariate and a constant the relationship between t and R2 is abs(t) = sqrt( R2/(1-R2) * (n-2) ), n = number of observations
#rVals <- seq(0.01,0.99,0.01)
#nObs <- 162
#tvals <- sqrt( rVals/(1-rVals) * (nObs-1) )
#plot( tvals ~ rVals ); abline(h=1.7, lwd=2,lty=2,col='red')
#1-pt(1.7,162)
#1-pt(tvals[10],162)

# to compare traveling wave analysis and modelling use the relationship between r (traveling wave analysis) and r2 (modelling)
# r = sqrt(r2)
  
  