rm( list=ls() );
mainDir <- '/analyse/Project0226/dataRepository'
targetDir <- '/analyse/Project0226/KastnerModel/data_KastnerClassic_pRF'
anatomiesDir_whole <- '/analyse/Project0226/KastnerModel/anatomies_KastnerClassic'
setwd( mainDir )

pRFDir <- dir( pattern='_Bars' )
pRFDir <- pRFDir[ c(1,2,3,5,7,8,10) ]
#setwd( 'anatomies' )
anatomiesDir <- dir( anatomiesDir_whole )
#anatomiesDir <- anatomiesDir[c(1,2,3,5,7,8,10)]

data.frame( pRFDir, anatomiesDir )

generalPurposeDir <- Sys.getenv( x='AFNI_TOOLBOXDIRGENERALPURPOSE' )
afniInstallDir <- Sys.getenv( x='AFNI_INSTALLDIR' ) 
source( sprintf('%s/scaleData.R', generalPurposeDir) )
source( sprintf('%s/AFNIio.R', afniInstallDir ) )
library( pracma )
library( abind )
library( circular )

for (nSubjects in 1:length(anatomiesDir)) { #subjId <- 1
  
  subjId <- nSubjects
  
  ############
  # pRF data #
  ############
  #setwd( mainDir )
  #setwd( kastnerDir[ subjId ] ); getwd()
  #system('rm -R zzz*')
  
  # clean data on the anatomy folder
  setwd( anatomiesDir_whole )
  setwd( anatomiesDir[ subjId ] ); getwd()
  system('rm -R ppp*')
  
  # back to the original epi folder (needs to be created)
  setwd( mainDir )
  setwd( pRFDir[ subjId ] ); getwd()
  
  #instr <- 'computeMeanTs.sh kastner_topUp/'; system( instr );
  #instr <- '3dTstat -prefix zzz_meanEpi4Coreg_kastner.nii.gz meanTs.nii'; system( instr )
  #system('rm meanTs.nii')
  
  #instr <- 'computeMeanTs.sh kastner_topUp_detrend/'; system( instr );
  #instr <- sprintf( 'mv meanTs.nii %s/%s/zzz_meanTs_kastner_topUp_dt.nii.gz', targetDir, kastnerDir[ subjId ] ); system( instr );
  instr <- sprintf('rm -R %s/%s', targetDir, pRFDir[ subjId ]); system( instr )
  instr <- sprintf('mkdir %s/%s', targetDir, pRFDir[ subjId ]); system( instr )
  instr <- sprintf( 'cp -R bars_topUp_detrend/ %s/%s/bars_topUp_detrend/', targetDir, pRFDir[ subjId ] ); system( instr );
  instr <- sprintf( 'cp -R coregistration_bars/ %s/%s/coregistration_bars/', targetDir, pRFDir[ subjId ] ); system( instr );
  instr <- sprintf( 'cp meanEpi4Coreg_bars.nii.gz %s/%s/meanEpi4Coreg_bars.nii.gz', targetDir, pRFDir[ subjId ] ); system( instr );
  
  # move to target dir
  setwd( targetDir  )
  setwd( pRFDir[ subjId ] ); getwd()
  
  #coregister
  setwd('coregistration_bars')
  system('sh coregistrationScript')
  setwd( '..' ); getwd()
  
  # gray matter mask
  instr <- sprintf( '3dAllineate -source %s/%s/volumetricData_layering_depth.nii.gz -prefix ppp_greyVolumetric_bars.nii.gz -1Dmatrix_apply coregistration_bars/invMat.1D -master meanEpi4Coreg_bars.nii.gz -final NN', anatomiesDir_whole, anatomiesDir[ subjId ] )
  system( instr )
  instr <- '3dcalc -a ppp_greyVolumetric_bars.nii.gz -expx \u0027step(a)*within(a,0.01,0.99)\u0027 -prefix ppp_grey_bars.nii.gz'; system( instr )

  # compute mean detrended time series
  instr <- 'computeMeanTs.sh bars_topUp_detrend/'; system( instr );
  instr <- 'mv meanTs.nii ppp_meanTs_bars_topUp_dt.nii.gz'; system( instr );
  
  # blur time series
  instr <- '3dBlurInMask -mask ppp_grey_bars.nii.gz -FWHM 1.5 -prefix ppp_meanTs_bars_topUp_dt_blur.nii.gz -input ppp_meanTs_bars_topUp_dt.nii.gz'; system( instr )
  
  # smooth time
  instr <- '3dTsmooth -prefix ppp_meanTs_bars_topUp_dt_blur_sm.nii.gz -lin ppp_meanTs_bars_topUp_dt_blur.nii.gz'; system( instr )
  
  # prf simple, prediction
  instr <- 'Rscript ../../code/generatePredictionsPrf_PSOCK_noMultParameter.R ppp_meanTs_bars_topUp_dt_blur_sm.nii.gz 4 0 0 0 ppp_savedPredictions_prfSimple'; system( instr )

  # prf simple, fit
  instr <- 'Rscript ../../code/fitPrf_PSOCK_noMultParameter.R ppp_meanTs_bars_topUp_dt_blur_sm.nii.gz ppp_grey_bars.nii.gz ppp_savedPredictions_prfSimple/ ppp_prf_nocss 0 0'; system( instr )
  
  # copy results in the anatomy, modelling
  instr <- sprintf( '3dAllineate -source ppp_prf_nocss_params.nii.gz -master ../../anatomies_KastnerClassic/%s/anatCopy.nii.gz -prefix ../../anatomies_KastnerClassic/%s/ppp_prf_nocss_params_anat.nii.gz -1Dmatrix_apply coregistration_bars/coregMat.1D -final NN', anatomiesDir[ subjId ], anatomiesDir[ subjId ] )
  system( instr )

  # move to anatomical directory and interpolate on surface
  setwd( sprintf('../../anatomies_KastnerClassic/%s/', anatomiesDir[ subjId ]) ); getwd()
  instr <- sprintf( '3dcalc -a ppp_prf_nocss_params_anat.nii.gz -exp \u0027a+0.0001\u0027 -prefix ppp_prf_nocss_params_anat_add.nii.gz'); system( instr )  
  instr <- sprintf( 'interpolateVolumeOnSurface.sh ppp_prf_nocss_params_anat_add.nii.gz surfaces_folder_left/'  ); system( instr )
  instr <- sprintf( 'interpolateVolumeOnSurface.sh ppp_prf_nocss_params_anat_add.nii.gz surfaces_folder_right/'  ); system( instr )
  instr <- sprintf( ' smoothSurfaceMap.sh ppp_prf_nocss_params_anat_add_interp_surfaces_folder_left/ 1.5 surfaces_folder_left/'  ); system( instr )
  instr <- sprintf( ' smoothSurfaceMap.sh ppp_prf_nocss_params_anat_add_interp_surfaces_folder_right/ 1.5 surfaces_folder_right/'  ); system( instr )
  
  #filesIn <- dir('ppp_prf_nocss_params_anat_add_interp_surfaces_folder_left_smoothed_1.5')
  #fileInTemp <- read.table( sprintf( 'ppp_prf_nocss_params_anat_add_interp_surfaces_folder_left_smoothed_1.5/%s', filesIn[1] ) )
  #polVals <- cart2pol( cbind( fileInTemp[,7], fileInTemp[,8] ) )
  #fileInTemp$newTheta <- polVals[,1]
  #fileInTemp$newEcc <- polVals[,2]
  
  #head( fileInTemp )
  #hist( fileInTemp[ fileInTemp[,18]>0.3 ,7] )
  #hist( polVals[ fileInTemp[,18]>0.3 ,2] )
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

