rm( list=ls() );
mainDir <- '/analyse/Project0226/KastnerModel'
setwd( mainDir )

dataDir <- 'data_KastnerModel'
dataFolders <- dir( dataDir )
anatomiesDir <- 'anatomies_KastnerModel'
anatomiesFolders <- dir( anatomiesDir )

dataFoldersAnalysis <- dataFolders[c(1,2,3,4,5,6,7,8)]
anatomiesFoldersAnalysis <- anatomiesFolders[c(1,2,3,4,5,6,7,8)]

data.frame( dataFoldersAnalysis, anatomiesFoldersAnalysis )# the name of the anatomy and the name of the data folder should correspond

generalPurposeDir <- Sys.getenv( x='AFNI_TOOLBOXDIRGENERALPURPOSE' )
afniInstallDir <- Sys.getenv( x='AFNI_INSTALLDIR' ) 
source( sprintf('%s/scaleData.R', generalPurposeDir) )
source( sprintf('%s/AFNIio.R', afniInstallDir ) )
library( pracma )
library( abind )
library( circular )

for ( subjId in 1:length(dataFoldersAnalysis) ) { # subjId <- 1; 
  
  ##########################
  # kastner modelling data #
  ##########################
  setwd( mainDir )
  setwd( dataDir ); setwd( dataFoldersAnalysis[ subjId ] ); getwd()
  system('rm -R zzz*')
  
  # move to anatomical directory and clean up
  setwd( sprintf('../../anatomies_KastnerModel/%s/', anatomiesFoldersAnalysis[ subjId ]) ); getwd()
  system('rm -R zzz_*')
  
  # back to epi directory
  setwd( mainDir )
  setwd( dataDir ); setwd( dataFoldersAnalysis[ subjId ] ); getwd()
  
  # gray matter mask
  instr <- sprintf( '3dAllineate -source ../../anatomies_KastnerModel/%s/volumetricData_layering_depth.nii.gz -prefix zzz_greyVolumetric_kastner.nii.gz -1Dmatrix_apply coreg/invMat.1D -master KASTMDL_stc_mc_tu_dt_mean.nii.gz -final NN', anatomiesFoldersAnalysis[ subjId ] )
  system( instr )
  instr <- '3dcalc -a zzz_greyVolumetric_kastner.nii.gz -expx \u0027step(a)*within(a,0.01,0.99)\u0027 -prefix zzz_grey_kastner.nii.gz'; system( instr )
  
  # modelling
  # blur in mask 
  instr <- '3dBlurInMask -mask zzz_grey_kastner.nii.gz -FWHM 0.001 -prefix zzz_KASTMDL_stc_mc_tu_dt_avg_blur.nii.gz -input KASTMDL_stc_mc_tu_dt_avg.nii.gz'; system( instr )

  # smooth time
  instr <- '3dTsmooth -prefix zzz_KASTMDL_stc_mc_tu_dt_avg_blur_smooth.nii.gz -lin zzz_KASTMDL_stc_mc_tu_dt_avg_blur.nii.gz'; system( instr )

  # kastner model prediction baseline model
  instr <- 'Rscript ../../code/generatePredictionsKastner_model_PSOCK.R zzz_KASTMDL_stc_mc_tu_dt_avg_blur_smooth.nii.gz zzz_kastner_model_baseline/ 0 baselineModel'; system( instr )

  # kastner model prediction kinematic model
  instr <- 'Rscript ../../code/generatePredictionsKastner_model_PSOCK.R zzz_KASTMDL_stc_mc_tu_dt_avg_blur_smooth.nii.gz zzz_kastner_model_kinematic/ 0 kinematic_baselineModel'; system( instr )
  
  # kastner model baseline, fit
  instr <- 'Rscript ../../code/fitEyeKinematics_PSOCK.R zzz_KASTMDL_stc_mc_tu_dt_avg_blur_smooth.nii.gz zzz_grey_kastner.nii.gz zzz_kastner_model_baseline/ zzz_kastner_model_baseline baselineModel'; system( instr )
  
  # kastner model kinematic, fit
  instr <- 'Rscript ../../code/fitEyeKinematics_PSOCK.R zzz_KASTMDL_stc_mc_tu_dt_avg_blur_smooth.nii.gz zzz_grey_kastner.nii.gz zzz_kastner_model_kinematic/ zzz_kastner_model_kinematic kinematic_baselineModel'; system( instr )
  
  # model difference AIC
  instr <- '3dcalc -a zzz_kastner_model_kinematic_params.nii.gz[12] -b zzz_kastner_model_baseline_params.nii.gz[12] -expr \u0027a-b\u0027 -prefix zzz_kastner_model_difference.nii.gz'; system( instr )
  
  # model difference adjR2
  instr <- '3dcalc -a zzz_kastner_model_kinematic_params.nii.gz[13] -b zzz_kastner_model_baseline_params.nii.gz[13] -expr \u0027a-b\u0027 -prefix zzz_kastner_model_difference_adjR2.nii.gz'; system( instr )
  
  #concatenate
  instr <- '3dTcat -prefix zzz_kastner_model_kinematic_params_cat.nii.gz zzz_kastner_model_kinematic_params.nii.gz zzz_kastner_model_difference.nii.gz zzz_kastner_model_difference_adjR2.nii.gz'; system( instr )
    
  # copy results in the anatomy
  instr <- sprintf( '3dAllineate -source zzz_kastner_model_kinematic_params_cat.nii.gz -master ../../anatomies_KastnerModel/%s/anatCopy.nii.gz -prefix ../../anatomies_KastnerModel/%s/zzz_kastner_model_params_anat.nii.gz -1Dmatrix_apply coreg/coregMat.1D -final NN', anatomiesFoldersAnalysis[ subjId ], anatomiesFoldersAnalysis[ subjId ] )
  system( instr )
  
  # move to anatomical directory and interpolate on surface
  setwd( sprintf('../../anatomies_KastnerModel/%s/', anatomiesFoldersAnalysis[ subjId ]) ); getwd()
  instr <- sprintf( '3dcalc -a zzz_kastner_model_params_anat.nii.gz -exp \u0027a+0.0001\u0027 -prefix zzz_kastner_model_params_anat_add.nii.gz'); system( instr )  
  instr <- sprintf( 'interpolateVolumeOnSurface.sh zzz_kastner_model_params_anat_add.nii.gz surfaces_folder_left/'  ); system( instr )
  instr <- sprintf( 'interpolateVolumeOnSurface.sh zzz_kastner_model_params_anat_add.nii.gz surfaces_folder_right/'  ); system( instr )
  
}  



#   
#   
#   
#     #kastner model, CCW
#   instr <- 'Rscript /analyse/Project0226/dataToSort/eyeMovPoster/codeFit/fitEyeKinematics_PSOCK.R zzz_CCW_mean_blur_sm.nii.gz zzz_mask_blur.nii.gz zzz_kastnerCCW/ zzz_CCW'; system( instr )
#   # #kastner model, CW, sinewave
#   # instr <- 'Rscript /analyse/Project0226/dataToSort/eyeMovPoster/codeFit/fitEyeKinematics_PSOCK.R zzz_CW_mean_blur_sm.nii.gz zzz_mask_blur.nii.gz zzz_kastner_sinewave/ zzz_CW_sinewave'; system( instr )
#   # #kastner model, CCW, sinewave
#   # instr <- 'Rscript /analyse/Project0226/dataToSort/eyeMovPoster/codeFit/fitEyeKinematics_PSOCK.R zzz_CCW_mean_blur_sm.nii.gz zzz_mask_blur.nii.gz zzz_kastner_sinewave/ zzz_CCW_sinewave'; system( instr )
#   #kastner model 2 components, CW both
#   instr <- 'Rscript /analyse/Project0226/dataToSort/eyeMovPoster/codeFit/fitEyeKinematics_PSOCK.R zzz_CW_mean_blur_sm.nii.gz zzz_mask_blur.nii.gz zzz_kastnerCW_2components_both/ zzz_CW_2components_both'; system( instr )
#   #kastner model 2 components, CCW both
#   instr <- 'Rscript /analyse/Project0226/dataToSort/eyeMovPoster/codeFit/fitEyeKinematics_PSOCK.R zzz_CCW_mean_blur_sm.nii.gz zzz_mask_blur.nii.gz zzz_kastnerCCW_2components_both/ zzz_CCW_2components_both'; system( instr )
# 
#   #kastner model 2 components, CW 10
#   instr <- 'Rscript /analyse/Project0226/dataToSort/eyeMovPoster/codeFit/fitEyeKinematics_PSOCK.R zzz_CW_mean_blur_sm.nii.gz zzz_mask_blur.nii.gz zzz_kastnerCW_2components_10/ zzz_CW_2components_10'; system( instr )
#   #kastner model 2 components, CCW 10
#   instr <- 'Rscript /analyse/Project0226/dataToSort/eyeMovPoster/codeFit/fitEyeKinematics_PSOCK.R zzz_CCW_mean_blur_sm.nii.gz zzz_mask_blur.nii.gz zzz_kastnerCCW_2components_10/ zzz_CCW_2components_10'; system( instr )
# 
#   #kastner model 2 components, CW 5
#   instr <- 'Rscript /analyse/Project0226/dataToSort/eyeMovPoster/codeFit/fitEyeKinematics_PSOCK.R zzz_CW_mean_blur_sm.nii.gz zzz_mask_blur.nii.gz zzz_kastnerCW_2components_05/ zzz_CW_2components_05'; system( instr )
#   #kastner model 2 components, CCW 5
#   instr <- 'Rscript /analyse/Project0226/dataToSort/eyeMovPoster/codeFit/fitEyeKinematics_PSOCK.R zzz_CCW_mean_blur_sm.nii.gz zzz_mask_blur.nii.gz zzz_kastnerCCW_2components_05/ zzz_CCW_2components_05'; system( instr )
#   
#   # instr <- 'compShiftPhase_model.sh zzz_CW_sinewave_params.nii.gz 33 zzz_CW_sinewave_params_shft33.nii.gz 6'; system( instr )
#   # instr <- 'compShiftPhase_model.sh zzz_CCW_sinewave_params.nii.gz 0 zzz_CCW_sinewave_params_shft33.nii.gz 6'; system( instr )
#   # instr <- 'compShiftPhase_model.sh zzz_CCW_sinewave_params_shft33.nii.gz 33 zzz_CCW_sinewave_params_shft33_inv.nii.gz 6'; system( instr )
#   
#   #set.seed(123); 
#   #par(mfrow=c(1,3)); 
#   #phv <- rnorm(10,0.5,0.9); 
#   #hist( phv, breaks=9, xlim=c(-3.14, 3.14) ); 
#   #phvt01 <- (phv+pi+0) ; 
#   #hist( phvt01, breaks=9, xlim=c(0, 6.28) );
#   #phvt02 <- (phv+0) %% (2*pi) ; 
#   #hist( phvt02, breaks=9, xlim=c(0, 6.28) );
#   
#   #kastner model, whole
#   #instr <- 'Rscript /analyse/Project0226/dataToSort/eyeMovPoster/codeFit/fitEyeKinematics.R zzz_whole_blur_smooth.nii.gz zzz_mask_blur.nii.gz zzz_kastner_whole/ zzz_whole_smooth'; system( instr )
#   
#   CWmod <- read.AFNI('zzz_CW_params.nii.gz')
#   CCWmod <- read.AFNI('zzz_CCW_params.nii.gz')
#   CWmodVol <- CWmod$brk
#   CCWmodVol <- CCWmod$brk
#   outVol <- array(0, dim( CCWmod$brk ) )
#   for (nVols in 1:dim( CCWmod$brk )[4] ) {
#     if (nVols!=6) { outVol[,,,nVols] <- ( CWmodVol[,,,nVols] + CCWmodVol[,,,nVols] ) / 2 }
#     if (nVols==6) { outVol[,,,nVols] <- atan2( ( sin( CWmodVol[,,,nVols] ) + sin( CCWmodVol[,,,nVols] ) ) / 2, 
#                                                ( cos( CWmodVol[,,,nVols] ) + cos( CCWmodVol[,,,nVols] ) ) / 2 )    }
#   }
#   write.AFNI( 'zzz_model_average.nii.gz', brk=outVol, origin=CWmod$origin, orient=CWmod$orient, defhead=CWmod$NI_head )
# 
#   CWmod <- read.AFNI('zzz_CW_2components_05_params.nii.gz')
#   CCWmod <- read.AFNI('zzz_CCW_2components_05_params.nii.gz')
#   CWmodVol <- CWmod$brk
#   CCWmodVol <- CCWmod$brk
#   outVol <- array(0, dim( CCWmod$brk ) )
#   for (nVols in 1:dim( CCWmod$brk )[4] ) {
#    if (nVols!=6) { outVol[,,,nVols] <- ( CWmodVol[,,,nVols] + CCWmodVol[,,,nVols] ) / 2 }
#    if (nVols==6) { outVol[,,,nVols] <- atan2( ( sin( CWmodVol[,,,nVols] ) + sin( CCWmodVol[,,,nVols] ) ) / 2, 
#                                               ( cos( CWmodVol[,,,nVols] ) + cos( CCWmodVol[,,,nVols] ) ) / 2 )    }
#   }
#   write.AFNI( 'zzz_model_average_2components_05.nii.gz', brk=outVol, origin=CWmod$origin, orient=CWmod$orient, defhead=CWmod$NI_head )
#   
#   CWmod <- read.AFNI('zzz_CW_2components_10_params.nii.gz')
#   CCWmod <- read.AFNI('zzz_CCW_2components_10_params.nii.gz')
#   CWmodVol <- CWmod$brk
#   CCWmodVol <- CCWmod$brk
#   outVol <- array(0, dim( CCWmod$brk ) )
#   for (nVols in 1:dim( CCWmod$brk )[4] ) {
#     if (nVols!=6) { outVol[,,,nVols] <- ( CWmodVol[,,,nVols] + CCWmodVol[,,,nVols] ) / 2 }
#     if (nVols==6) { outVol[,,,nVols] <- atan2( ( sin( CWmodVol[,,,nVols] ) + sin( CCWmodVol[,,,nVols] ) ) / 2, 
#                                                ( cos( CWmodVol[,,,nVols] ) + cos( CCWmodVol[,,,nVols] ) ) / 2 )    }
#   }
#   write.AFNI( 'zzz_model_average_2components_10.nii.gz', brk=outVol, origin=CWmod$origin, orient=CWmod$orient, defhead=CWmod$NI_head )
#   
#   CWmod <- read.AFNI('zzz_CW_2components_both_params.nii.gz')
#   CCWmod <- read.AFNI('zzz_CCW_2components_both_params.nii.gz')
#   CWmodVol <- CWmod$brk
#   CCWmodVol <- CCWmod$brk
#   outVol <- array(0, dim( CCWmod$brk ) )
#   for (nVols in 1:dim( CCWmod$brk )[4] ) {
#     if (nVols!=6) { outVol[,,,nVols] <- ( CWmodVol[,,,nVols] + CCWmodVol[,,,nVols] ) / 2 }
#     if (nVols==6) { outVol[,,,nVols] <- atan2( ( sin( CWmodVol[,,,nVols] ) + sin( CCWmodVol[,,,nVols] ) ) / 2, 
#                                                ( cos( CWmodVol[,,,nVols] ) + cos( CCWmodVol[,,,nVols] ) ) / 2 )    }
#   }
#   write.AFNI( 'zzz_model_average_2components_both.nii.gz', brk=outVol, origin=CWmod$origin, orient=CWmod$orient, defhead=CWmod$NI_head )
#   
#   
# 
# }
# 
# # simple model with one covariate and a constant the relationship between t and R2 is abs(t) = sqrt( R2/(1-R2) * (n-2) ), n = number of observations
# #rVals <- seq(0.01,0.99,0.01)
# #nObs <- 162
# #tvals <- sqrt( rVals/(1-rVals) * (nObs-1) )
# #plot( tvals ~ rVals ); abline(h=1.7, lwd=2,lty=2,col='red')
# #1-pt(1.7,162)
# #1-pt(tvals[10],162)
# 
# # to compare traveling wave analysis and modelling use the relationship between r (traveling wave analysis) and r2 (modelling)
# # r = sqrt(r2)
#   
#   
