rm( list=ls() );
mainDir <- '/analyse/Project0226/dataRepository'
targetDir_pRF <- '/analyse/Project0226/KastnerModel/data_KastnerClassic_pRF'
targetDir_KastnerClassic <- '/analyse/Project0226/KastnerModel/data_KastnerClassic'
anatomiesDir_whole <- '/analyse/Project0226/KastnerModel/anatomies_KastnerClassic'
resultsDir <- '/analyse/Project0226/KastnerModel/results_KastnerClassic'
filenameOut <- 'testResults_all_in_Kastner_classic_withHrf_08102021.RData'
setwd( mainDir )

pRFDir <- dir( targetDir_pRF )
KCDir <- dir( targetDir_KastnerClassic )
anatomiesDir <- dir( anatomiesDir_whole )

data.frame( pRFDir, KCDir, anatomiesDir )

generalPurposeDir <- Sys.getenv( x='AFNI_TOOLBOXDIRGENERALPURPOSE' )
afniInstallDir <- Sys.getenv( x='AFNI_INSTALLDIR' ) 
source( sprintf('%s/scaleData.R', generalPurposeDir) )
source( sprintf('%s/AFNIio.R', afniInstallDir ) )
library( pracma )
library( abind )
library( circular )

# bring pRF data into kastner classic data
setwd( anatomiesDir_whole ); dir()
for ( iAnat in 1:length( anatomiesDir ) ) { #iAnat <- 1
  setwd( anatomiesDir_whole )
  setwd( anatomiesDir[ iAnat ] )
  if ( file.exists( sprintf( '../../data_KastnerClassic/%s/zzz_params_prf.nii.gz', KCDir[ iAnat ] ) ) ) { system( sprintf( 'rm ../../data_KastnerClassic/%s/zzz_params_prf.nii.gz', KCDir[ iAnat ] ) ) }
  instr <- sprintf( '3dAllineate -source ppp_prf_nocss_params_anat_add.nii.gz -prefix ../../data_KastnerClassic/%s/zzz_params_prf.nii.gz -1Dmatrix_apply ../../data_KastnerClassic/%s/coregistration_kastner/invMat.1D -master ../../data_KastnerClassic/%s/meanEpi4Coreg_kastner.nii.gz -final NN', KCDir[ iAnat ], KCDir[ iAnat ], KCDir[ iAnat ] )
  system( instr )
}


###################################################################################
# extract model data
extractModelData <- function( targetDir, subjDir, modelFileName, modelDFName ) { #targetDir <- targetDir_KastnerClassic; subjDir <- KCDir; modelFileName <- 'zzz_params_prf.nii.gz'; modelDFName <- 'pRF'
  flagOut <- 1
  for (nSubjects in 1:length(anatomiesDir)) { #nSubjects <- 1
    
    subjId <- nSubjects
    
    # set the Kastner classic folder
    setwd( targetDir ); getwd()
    setwd( subjDir[ subjId ] ); getwd()
    modelParams <- read.AFNI( modelFileName )
    modelParamsVolume <- modelParams$brk
    
    # get data from rois
    setwd('roi_epi/')
    roiFiles <- dir( pattern='*.BRIK' )
    print( roiFiles )
    print('..........................')
    print('..........................')
    print('rois to process...........')
    print('..........................')
    print('..........................')
    
    for ( nRois in 1:length(roiFiles)) { #length(roiFiles) nRois <- 1
      
      print('.............................')
      print('.............................')
      print('.............................')
      print( sprintf('processing roi: %s, participant: %s; ', roiFiles[nRois], subjDir[ subjId ]) )
      print('.............................')
      print('.............................')
      print('.............................')
      
      roiTemp <- read.AFNI( roiFiles[ nRois ] )
      roiTempVolume <- roiTemp$brk
      roiTempIdx <- which( roiTempVolume==1 ) 
      df_roi <- data.frame( rep( modelDFName, length( roiTempIdx ) ), rep( subjDir[ subjId ], length( roiTempIdx ) ), rep( roiFiles[ nRois ], length( roiTempIdx ) ) , array(0, c( length( roiTempIdx ), dim(modelParamsVolume)[4]  ) ) )
      
      for ( nElements in 1:dim(modelParamsVolume)[4] ) { #pRFElements <- 1
        currentModelData <- modelParamsVolume[,,,nElements]
        df_roi[, (nElements+3) ] <- currentModelData[ roiTempIdx ]
      }
      
      if (flagOut==1) { df_roi_out <- df_roi; flagOut <- 0 }
      if (flagOut==0) { df_roi_out <- rbind( df_roi_out, df_roi )  }
      
    }
    
  }
  return( df_roi_out )
}
#df_prf <- extractModelData( targetDir_pRF, pRFDir, 'ppp_prf_nocss_params.nii.gz', 'pRF' )

df_prf <- extractModelData( targetDir_KastnerClassic, KCDir, 'zzz_params_prf.nii.gz', 'pRF' )

df_modelBoth <- extractModelData( targetDir_KastnerClassic, KCDir, 'zzz_model_wAverage_both.nii.gz', 'model_both' )
df_modelBoth_CW <- extractModelData( targetDir_KastnerClassic, KCDir, 'zzz_CW_2components_both_params.nii.gz', 'model_both_CW' )
df_modelBoth_CCW <- extractModelData( targetDir_KastnerClassic, KCDir, 'zzz_CCW_2components_both_params.nii.gz', 'model_both_CCW' )

df_model05 <- extractModelData( targetDir_KastnerClassic, KCDir, 'zzz_model_wAverage_05.nii.gz', 'model_05' )
df_model05_CW <- extractModelData( targetDir_KastnerClassic, KCDir, 'zzz_CW_2components_05_params.nii.gz', 'model_05_CW' )
df_model05_CCW <- extractModelData( targetDir_KastnerClassic, KCDir, 'zzz_CCW_2components_05_params.nii.gz', 'model_05_CCW' )

df_model10 <- extractModelData( targetDir_KastnerClassic, KCDir, 'zzz_model_wAverage_10.nii.gz', 'model_10' )
df_model10_CW <- extractModelData( targetDir_KastnerClassic, KCDir, 'zzz_CW_2components_10_params.nii.gz', 'model_10_CW' )
df_model10_CCW <- extractModelData( targetDir_KastnerClassic, KCDir, 'zzz_CCW_2components_10_params.nii.gz', 'model_10_CCW' )

df_phase05 <- extractModelData( targetDir_KastnerClassic, KCDir, 'zzz_phase_encoded_wAverage_05.nii.gz', 'phase_05' )
df_phase05_CW <- extractModelData( targetDir_KastnerClassic, KCDir, 'zzz_trWave05_CW.nii.gz', 'phase_05_CW' )
df_phase05_CCW <- extractModelData( targetDir_KastnerClassic, KCDir, 'zzz_trWave05_CCW.nii.gz', 'phase_05_CCW' )

df_phase10 <- extractModelData( targetDir_KastnerClassic, KCDir, 'zzz_phase_encoded_wAverage_10.nii.gz', 'phase_10' )
df_phase10_CW <- extractModelData( targetDir_KastnerClassic, KCDir, 'zzz_trWave10_CW.nii.gz', 'phase_10_CW' )
df_phase10_CCW <- extractModelData( targetDir_KastnerClassic, KCDir, 'zzz_trWave10_CCW.nii.gz', 'phase_10_CCW' )

names( df_prf ) <- c('exp','subj','roiName','xPos','yPos','sigmaPos','sigmaNeg','hrf_a1','hrf_a2','multPar','staticNonLin','nCycles','angle','phase','varExp','intercept','slope','logLik','AIC','theta','radius','fwhmCenter','surroundSize')

names( df_modelBoth ) <- c('exp','subj','roiName','hrf_a1','hrf_a2','hrf_b1','hrf_b2','hrf_c','mu','kappa','nonLin','varExp','intercept','slope','logLik','AIC','adjR2')
names( df_modelBoth_CW ) <- c('exp','subj','roiName','hrf_a1','hrf_a2','hrf_b1','hrf_b2','hrf_c','mu','kappa','nonLin','varExp','intercept','slope','logLik','AIC','adjR2')
names( df_modelBoth_CCW ) <- c('exp','subj','roiName','hrf_a1','hrf_a2','hrf_b1','hrf_b2','hrf_c','mu','kappa','nonLin','varExp','intercept','slope','logLik','AIC','adjR2')

names( df_model05 ) <- c('exp','subj','roiName','hrf_a1','hrf_a2','hrf_b1','hrf_b2','hrf_c','mu','kappa','nonLin','varExp','intercept','slope','logLik','AIC','adjR2')
names( df_model05_CW ) <- c('exp','subj','roiName','hrf_a1','hrf_a2','hrf_b1','hrf_b2','hrf_c','mu','kappa','nonLin','varExp','intercept','slope','logLik','AIC','adjR2')
names( df_model05_CCW ) <- c('exp','subj','roiName','hrf_a1','hrf_a2','hrf_b1','hrf_b2','hrf_c','mu','kappa','nonLin','varExp','intercept','slope','logLik','AIC','adjR2')

names( df_model10 ) <- c('exp','subj','roiName','hrf_a1','hrf_a2','hrf_b1','hrf_b2','hrf_c','mu','kappa','nonLin','varExp','intercept','slope','logLik','AIC','adjR2')
names( df_model10_CW ) <- c('exp','subj','roiName','hrf_a1','hrf_a2','hrf_b1','hrf_b2','hrf_c','mu','kappa','nonLin','varExp','intercept','slope','logLik','AIC','adjR2')
names( df_model10_CCW ) <- c('exp','subj','roiName','hrf_a1','hrf_a2','hrf_b1','hrf_b2','hrf_c','mu','kappa','nonLin','varExp','intercept','slope','logLik','AIC','adjR2')

names( df_phase05 ) <- c('exp','subj','roiName','amp','co','phase')
names( df_phase05_CW ) <- c('exp','subj','roiName','amp','co','phase')
names( df_phase05_CCW ) <- c('exp','subj','roiName','amp','co','phase')

names( df_phase10 ) <- c('exp','subj','roiName','amp','co','phase')
names( df_phase10_CW ) <- c('exp','subj','roiName','amp','co','phase')
names( df_phase10_CCW ) <- c('exp','subj','roiName','amp','co','phase')

save( df_prf, df_modelBoth, df_modelBoth_CW, df_modelBoth_CCW, 
      df_model05, df_model05_CW, df_model05_CCW, 
      df_model10, df_model10_CW, df_model10_CCW, 
      df_phase05, df_phase05_CW, df_phase05_CCW, 
      df_phase10, df_phase10_CW, df_phase10_CCW,
      file=sprintf('%s/%s', resultsDir, filenameOut )  )

