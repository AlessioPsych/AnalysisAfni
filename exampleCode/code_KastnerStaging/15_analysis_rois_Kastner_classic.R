rm( list=ls() );
mainDir <- '/analyse/Project0226/KastnerModel/staging_area_Kastner'
KastnerDir <- 'KastnerClassic_hrf_modulation'
pRFDir <- 'pRF'
figuresOutput <- 'figures_KastnerClassic_hrf_modulation'
resultsDir <- 'results'
setwd( mainDir )
setwd( resultsDir )

generalPurposeDir <- Sys.getenv( x='AFNI_TOOLBOXDIRGENERALPURPOSE' )
afniInstallDir <- Sys.getenv( x='AFNI_INSTALLDIR' ) 
source( sprintf('%s/scaleData.R', generalPurposeDir) )
source( sprintf('%s/AFNIio.R', afniInstallDir ) )
library( pracma )
library( abind )
library( circular )
library( lme4 )
library( lmerTest )
library( ggplot2 )
library( pROC )

#### current folder ####
getwd()

#### load data, phase 05 ####
load( sprintf('%s/_phase_encoded_wAverage_05.nii.gz.RData', KastnerDir ) )
df_phase_encoded_wAverage_05 <- outModelDf
rm( list = c( 'outModelDf' ) )
str( df_phase_encoded_wAverage_05 )

#### load data, phase 05 CW ####
load( sprintf('%s/_trWave05_CW.nii.gz.RData', KastnerDir ) )
df_phase_encoded_CW_05 <- outModelDf
rm( list = c( 'outModelDf' ) )
str( df_phase_encoded_CW_05 )

#### load data, phase 05 CCW ####
load( sprintf('%s/_trWave05_CCW.nii.gz.RData', KastnerDir ) )
df_phase_encoded_CCW_05 <- outModelDf
rm( list = c( 'outModelDf' ) )
str( df_phase_encoded_CCW_05 )




#### load data, phase 10 ####
load( sprintf('%s/_phase_encoded_wAverage_10.nii.gz.RData', KastnerDir ) )
df_phase_encoded_wAverage_10 <- outModelDf
rm( list = c( 'outModelDf' ) )
str( df_phase_encoded_wAverage_10 )

#### load data, phase 10 CW ####
load( sprintf('%s/_trWave10_CW.nii.gz.RData', KastnerDir ) )
df_phase_encoded_CW_10 <- outModelDf
rm( list = c( 'outModelDf' ) )
str( df_phase_encoded_CW_10 )

#### load data, phase 10 CCW ####
load( sprintf('%s/_trWave10_CCW.nii.gz.RData', KastnerDir ) )
df_phase_encoded_CCW_10 <- outModelDf
rm( list = c( 'outModelDf' ) )
str( df_phase_encoded_CCW_10 )





#### load data, model 05 ####
load( sprintf('%s/_model_wAverage_05.nii.gz.RData', KastnerDir ) )
df_model_wAverage_05 <- outModelDf
rm( list = c( 'outModelDf' ) )
str( df_model_wAverage_05 )

#### load data, phase 05 CW ####
load( sprintf('%s/_CW_2components_05_params.nii.gz.RData', KastnerDir ) )
df_model_05_CW <- outModelDf
rm( list = c( 'outModelDf' ) )
str( df_model_05_CW )

#### load data, phase 05 CCW ####
load( sprintf('%s/_CCW_2components_05_params.nii.gz.RData', KastnerDir ) )
df_model_05_CCW <- outModelDf
rm( list = c( 'outModelDf' ) )
str( df_model_05_CCW )




#### load data, model 10 ####
load( sprintf('%s/_model_wAverage_10.nii.gz.RData', KastnerDir ) )
df_model_wAverage_10 <- outModelDf
rm( list = c( 'outModelDf' ) )
str( df_model_wAverage_10 )

#### load data, phase 10 CW ####
load( sprintf('%s/_CW_2components_10_params.nii.gz.RData', KastnerDir ) )
df_model_10_CW <- outModelDf
rm( list = c( 'outModelDf' ) )
str( df_model_10_CW )

#### load data, phase 10 CCW ####
load( sprintf('%s/_CCW_2components_10_params.nii.gz.RData', KastnerDir ) )
df_model_10_CCW <- outModelDf
rm( list = c( 'outModelDf' ) )
str( df_model_10_CCW )




#### load data, model both ####
load( sprintf('%s/_model_wAverage_both.nii.gz.RData', KastnerDir ) )
df_model_wAverage_both <- outModelDf
rm( list = c( 'outModelDf' ) )
str( df_model_wAverage_both )

#### load data, phase both CW ####
load( sprintf('%s/_CW_2components_both_params.nii.gz.RData', KastnerDir ) )
df_model_both_CW <- outModelDf
rm( list = c( 'outModelDf' ) )
str( df_model_both_CW )

#### load data, phase both CCW ####
load( sprintf('%s/_CCW_2components_both_params.nii.gz.RData', KastnerDir ) )
df_model_both_CCW <- outModelDf
rm( list = c( 'outModelDf' ) )
str( df_model_both_CCW )




#### load data, model pRF ####
load( sprintf('%s/_pRF_params.nii.gz.RData', pRFDir ) )
df_model_pRF <- outModelDf
rm( list = c( 'outModelDf' ) )
str( df_model_pRF )

table( df_model_pRF$participant, df_model_pRF$ROI_name )
table( df_model_wAverage_05$participant, df_model_wAverage_05$ROI_name )



# adjusted r squared for prf data
#nParams_prf <- 3
#df_prf_filt$adjR2 <- 1 - (1-df_prf_filt$varExp) * ( (155-1) / (155-nParams_prf-1) )

df_model_pRF$radius
df_model_pRF$participant
df_model_pRF$ROI_name_basic
df_model_pRF$varExp

# cut data function
dataIn <- df_model_pRF
participantId <- 'AFH28'
propCut <- c(0,1,0.1)
roiInd <- 'V1'
varThr <- 0.4

cutData <- function( dataIn, participantId, propCut, roiInd, varThr ) {
  dataSel <- dataIn[ dataIn$participant==participantId & 
                       dataIn$ROI_name_basic==roiInd, ]
  varCut <- dataSel$varExp > quantile( dataSel$varExp, varThr )
  
  keptIdx <- which( varCut )
  dataSel <- dataSel[ varCut, ]
  eccCut <- .bincode( dataSel$radius, 
                      breaks = quantile( dataSel$radius, probs =  seq( propCut[1], propCut[2], propCut[3] ) ), 
                      include.lowest = TRUE )
  xCut <- .bincode( dataSel$`x-Pos`, 
                      breaks = quantile( dataSel$`x-Pos`, probs =  seq( propCut[1], propCut[2], propCut[3] ) ), 
                      include.lowest = TRUE )
  yCut <- .bincode( dataSel$`y-Pos`, 
                    breaks = quantile( dataSel$`y-Pos`, probs =  seq( propCut[1], propCut[2], propCut[3] ) ), 
                    include.lowest = TRUE )
  dataSel$eccCut <- as.numeric( eccCut )
  dataSel$xCut <- as.numeric( xCut )
  dataSel$yCut <- as.numeric( yCut )
  dataSel$keptIdx <- keptIdx
  return( dataSel )
}

##### for each participant and roi, cut the data #####
flagDf <- 1
for ( participantId in unique( df_model_pRF$participant ) ) {
  df_prf_filt_subj <- df_model_pRF[ df_model_pRF$participant==participantId, ]
  for ( roiId in unique( df_model_pRF$ROI_name_basic ) ) {
    dataDfTemp <- cutData(df_model_pRF, participantId, c(0,1,0.05), roiId, 0.6  ) 
    if (flagDf==1) { df_prf_filt_clean_full <- dataDfTemp; flagDf <- 0  }
    if (flagDf==0) { df_prf_filt_clean_full <- rbind( df_prf_filt_clean_full, dataDfTemp ) }
  }
}
str( df_prf_filt_clean_full ) 
table( df_prf_filt_clean_full$ROI_name_basic, df_prf_filt_clean_full$participant )

##### select only ROIs of interest #####

df_prf_filt_clean_full$ROI_name_basic <- as.character( df_prf_filt_clean_full$ROI_name_basic )
df_prf_filt_clean <- subset( df_prf_filt_clean_full, ROI_name_basic=='V1' | ROI_name_basic=='V2' | ROI_name_basic=='V3'  | ROI_name_basic=='hV4' 
                                       | ROI_name_basic=='MST' | ROI_name_basic=='hMT' | ROI_name_basic=='V3a' | ROI_name_basic=='V3b' | ROI_name_basic=='IPS0' | ROI_name_basic=='IPS1'
                                       | ROI_name_basic=='IPS2' | ROI_name_basic=='IPS3' | ROI_name_basic=='IPS4' )
df_prf_filt_clean$ROI_name_basic <- ifelse( df_prf_filt_clean$ROI_name_basic == 'hMT' | df_prf_filt_clean$ROI_name_basic == 'MST', 'hMT/MST', df_prf_filt_clean$ROI_name_basic )
roiLevels <- unique(  df_prf_filt_clean$ROI_name_basic )
df_prf_filt_clean$roi_sorted <- factor( df_prf_filt_clean$ROI_name_basic, levels=roiLevels ) 
df_prf_filt_clean$ROI_name_basic <- df_prf_filt_clean$roi_sorted

##### size x eccentricy, single subjects, plot and analysis ##### 

# size x eccentricy, single subjects
size_ecc_df_base <- aggregate( sigmaPos ~ eccCut * roi_sorted * participant, mean, data=df_prf_filt_clean )
size_ecc_df_base$radius <- aggregate( radius ~ eccCut * roi_sorted * participant, mean, data=df_prf_filt_clean )$radius
size_ecc_df_base$varExp <- aggregate( varExp ~ eccCut * roi_sorted * participant, mean, data=df_prf_filt_clean )$varExp

# size x eccentricity, global average
size_ecc_df_base_globalAverage <- aggregate( sigmaPos ~ eccCut * roi_sorted, mean, data=size_ecc_df_base) 
size_ecc_df_base_globalAverage$radius <- aggregate( radius ~ eccCut * roi_sorted, mean, data=size_ecc_df_base)$radius

# generate figure
x11( width=9, height=9 )

par(mfrow=c(4,5))
subjCols <- rainbow( length( unique( df_prf_filt_clean$roi_sorted ) ), star=0, end=0.8 )
for ( participantId in unique( df_prf_filt_clean$participant ) ) {
  plot( sigmaPos ~ radius, data=size_ecc_df_base[ size_ecc_df_base$roi_sorted=='V1' & size_ecc_df_base$participant==participantId, ], 
        col=subjCols[ 1 ], xlab='eccentricity (dva)', ylab='size (dva)', cex.lab=1.55, cex.axis=1.5,
        ylim=c(0,10), xlim=c(0,10), pch=8, las=1, bty='n' ); 
  abline( lm( sigmaPos ~ radius, data=size_ecc_df_base[ size_ecc_df_base$roi_sorted=='V1' & size_ecc_df_base$participant==participantId, ] ), col=subjCols[ 1 ] )
  text( 2.5, 9, participantId, cex=1.25 )
  roiCounter <- 2
  for (roiIdx in roiLevels ) {
    dataPlot <- size_ecc_df_base[ size_ecc_df_base$roi_sorted==roiIdx & size_ecc_df_base$participant==participantId, ]; points( dataPlot$radius, dataPlot$sigmaPos, col=subjCols[ roiCounter ], pch=8 ); abline( lm( dataPlot$sigmaPos ~ dataPlot$radius ), col=subjCols[ roiCounter ] )
    roiCounter <- roiCounter + 1
  }
}

# size x eccentricity, global average, plot
plot( sigmaPos ~ radius, data=size_ecc_df_base_globalAverage[ size_ecc_df_base_globalAverage$roi_sorted=='V2', ], col=subjCols[ 1 ], 
      ylim=c(0,10), xlim=c(0,10), xlab='eccentricity (dva)', ylab='size (dva)', cex.lab=1.55, cex.axis=1.5,
      pch=8, las=1, bty='n' ); abline( lm( sigmaPos ~ radius, data=size_ecc_df_base_globalAverage[ size_ecc_df_base_globalAverage$roi_sorted=='V2', ] ), col=subjCols[ 1 ] )
text( 3.2, 9, 'Average', cex=1.5 )
roiCounter <- 2
for (roiIdx in roiLevels ) {
  dataPlot <- size_ecc_df_base_globalAverage[ size_ecc_df_base_globalAverage$roi_sorted==roiIdx, ]; points( dataPlot$radius, dataPlot$sigmaPos, col=subjCols[ roiCounter ], pch=8 ); abline( lm( dataPlot$sigmaPos ~ dataPlot$radius ), col=subjCols[ roiCounter ] )
  roiCounter <- roiCounter + 1
}

# legend plot
plot(0,axes=FALSE,bty='n',xlab='',ylab='',type='n',xlim=c(0,1),ylim=c(0,1))
legend(0, 1.0,legend=roiLevels[1:5],fill=subjCols[1:5],bty='n',cex=0.8)
legend(0.4, 1.0,legend=roiLevels[5:9],fill=subjCols[5:9],bty='n',cex=0.8)
legend(0.8, 1.0,legend=roiLevels[10:12],fill=subjCols[10:12],bty='n',cex=0.8)

dev.copy2pdf( file=sprintf('%s/%s/prf_sizeEccentricity.pdf', mainDir, figuresOutput ), width=9, height=9 )
graphics.off()

#### statistical analysis ####
size_ecc_df_base$roi_factor <- as.factor( size_ecc_df_base$roi_sorted )
size_ecc_df_base$eccCut_scalar <- as.numeric( size_ecc_df_base$eccCut )
size_ecc_df_base$subject_factor <- as.factor( size_ecc_df_base$participant )
summary( lmer( sigmaPos ~ eccCut_scalar * roi_factor + ( 1 | subject_factor ), data=size_ecc_df_base ) )

##### prf variance explained plot and analysis, adj.r2 ##### 
x11(width=4.5, height=4.5)

par(mfrow=c(1,1))
size_ecc_df_base_globalAverage_var <- aggregate( varExp ~ participant * roi_sorted, mean, data=df_prf_filt_clean )  #here varExp refers to adjR2, see line 88 and line 90
bp <- boxplot( size_ecc_df_base_globalAverage_var$varExp ~ size_ecc_df_base_globalAverage_var$roi_sorted, ylim=c(0,1), 
         frame=FALSE, las=2, col=subjCols, outline=FALSE, xlab='ROI', ylab='adj. R^2 (proportion)', cex.lab=1.35, cex.axis=1.25 )
sMod <- summary( lmer( varExp ~ roi_sorted + ( 1 | participant ), data=size_ecc_df_base_globalAverage_var ) )
for (i in 1:length(sMod$coefficients[,5])) { 
  if ( sMod$coefficients[i,5] < 0.001 ) {
    if (i==1) { text( i, bp$stats[5,i]+0.1, '*', cex=3, col='red' ) }
    if (i!=1) { text( i, bp$stats[5,i]+0.1, '*', cex=3 ) }
  }
}

dev.copy2pdf( file=sprintf('%s/%s/prf_varExp_adjR2.pdf', mainDir, figuresOutput ), width=4.5, height=4.5 )
graphics.off()

##### prf variance explained plot and analysis, adj.r2, quantiles ##### 
subsetDF_quantile <- c()
for ( p in unique( df_prf_filt_clean$participant ) ) { # p <- 'AFH28'
    subsetDF <- subset( df_prf_filt_clean, participant==p )
    subsetDF$varExp_quantile <- as.numeric( cut( subsetDF$varExp, 
                                             breaks=quantile( subsetDF$varExp, seq( 0, 1, 0.01 ) ), 
                                             include.lowest = TRUE ) ) / 100
    subsetDF_quantile <- rbind( subsetDF_quantile, subsetDF )
}

x11(width=4.5, height=4.5)

par(mfrow=c(1,1))
size_ecc_df_base_globalAverage_var <- aggregate( varExp_quantile ~ participant * roi_sorted, mean, data=subsetDF_quantile )  
bp <- boxplot( size_ecc_df_base_globalAverage_var$varExp ~ size_ecc_df_base_globalAverage_var$roi_sorted, ylim=c(0,1), 
               frame=FALSE, las=2, col=subjCols, outline=FALSE, xlab='ROI', ylab='adj. R^2 (proportion)', cex.lab=1.35, cex.axis=1.25 )
sMod <- summary( lmer( varExp_quantile ~ roi_sorted + ( 1 | participant ), data=size_ecc_df_base_globalAverage_var ) )
for (i in 1:length(sMod$coefficients[,5])) { 
  if ( sMod$coefficients[i,5] < 0.001 ) {
    if (i==1) { text( i, bp$stats[5,i]+0.1, '*', cex=3, col='red' ) }
    if (i!=1) { text( i, bp$stats[5,i]+0.1, '*', cex=3 ) }
  }
}

dev.copy2pdf( file=sprintf('%s/%s/prf_varExp_adjR2_quantiles.pdf', mainDir, figuresOutput ), width=4.5, height=4.5 )
graphics.off()

##### model kastner (filtered on model_kastner classic data itself), prepare data ##### 
str(df_model_wAverage_05)

#df_phase_encoded_wAverage_05
#df_model_wAverage_05
df_kastner_model <- df_model_wAverage_05; 

df_kastner_model$logLik_05 <- df_kastner_model$logLik; nParams_05 <- 4 # location / width / intercept / slope
df_kastner_model$AIC_05 <- 2*nParams_05 - 2*df_kastner_model$logLik_05
df_kastner_model$adjR2_05 <- 1 - (1-df_kastner_model$varExp) * ( (279-1) / (279-nParams_05-1) )
df_kastner_model$mu05_CW <- df_model_05_CW$mu #df_model05_CW$mu
df_kastner_model$mu05_CCW <- df_model_05_CCW$mu
df_kastner_model$kappa05_CW <- df_model_05_CW$kappa
df_kastner_model$kappa05_CCW <- df_model_05_CCW$kappa
df_kastner_model$adjR2_05_CW <- 1 - (1-df_model_05_CW$varExp) * ( (279-1) / (279-nParams_05-1) )
df_kastner_model$adjR2_05_CCW <- 1 - (1-df_model_05_CCW$varExp) * ( (279-1) / (279-nParams_05-1) )
df_kastner_model$adjR2_05_CW_CCW_avg <- ( df_kastner_model$adjR2_05_CW + df_kastner_model$adjR2_05_CCW ) / 2

df_kastner_model$logLik_10 <- df_model_wAverage_10$logLik; nParams_10 <- 4 # location / width / intercept / slope
df_kastner_model$AIC_10 <- 2*nParams_10 - 2*df_kastner_model$logLik_10 
df_kastner_model$adjR2_10 <- 1 - (1-df_model_wAverage_10$varExp) * ( (279-1) / (279-nParams_10-1) )
df_kastner_model$mu10_CW <- df_model_10_CW$mu
df_kastner_model$mu10_CCW <- df_model_10_CCW$mu
df_kastner_model$kappa10_CW <- df_model_10_CW$kappa
df_kastner_model$kappa10_CCW <- df_model_10_CCW$kappa
df_kastner_model$adjR2_10_CW <- 1 - (1-df_model_10_CW$varExp) * ( (279-1) / (279-nParams_10-1) )
df_kastner_model$adjR2_10_CCW <- 1 - (1-df_model_10_CCW$varExp) * ( (279-1) / (279-nParams_10-1) )

df_kastner_model$logLik_both <- df_model_wAverage_both$logLik; nParams_both <- 5 # location / width / intercept / slope / ratio
df_kastner_model$AIC_both <- 2*nParams_both - 2*df_kastner_model$logLik_both
df_kastner_model$adjR2_both <- 1 - (1-df_model_wAverage_both$varExp) * ( (279-1) / (279-nParams_both-1) )
df_kastner_model$both_ratio <- df_model_wAverage_both$nLin
df_kastner_model$muBoth_CW <- df_model_both_CW$mu
df_kastner_model$muBoth_CCW <- df_model_both_CCW$mu
df_kastner_model$kappaBoth_CW <- df_model_both_CW$kappa
df_kastner_model$kappaBoth_CCW <- df_model_both_CCW$kappa
df_kastner_model$adjR2_both_CW <- 1 - (1-df_model_both_CW$varExp) * ( (279-1) / (279-nParams_both-1) )
df_kastner_model$adjR2_both_CCW <- 1 - (1-df_model_both_CCW$varExp) * ( (279-1) / (279-nParams_both-1) )

df_phase05 <- df_phase_encoded_wAverage_05
df_phase10 <- df_phase_encoded_wAverage_10
df_phase05_CW <- df_phase_encoded_CW_05
df_phase05_CCW <- df_phase_encoded_CCW_05
df_phase10_CW <- df_phase_encoded_CW_10
df_phase10_CCW <- df_phase_encoded_CCW_10
nParams_phase05 <- 2 # phase / amplitude 
nParams_phase10 <- 2 # phase / amplitude 
df_kastner_model$adjR2_phase05 <- 1 - ( 1-df_phase05$co^2 ) * ( (279-1) / (279-nParams_phase05-1) )
df_kastner_model$adjR2_phase10 <- 1 - ( 1-df_phase10$co^2 ) * ( (279-1) / (279-nParams_phase10-1) )
df_kastner_model$phase05 <- df_phase05$phase
df_kastner_model$phase05_CW <- ( df_phase05_CW$ph - 0.62 ) %% ( 2*pi )
df_kastner_model$phase05_CCW <- ( df_phase05_CCW$ph - 0.62 ) %% ( 2*pi )
df_kastner_model$phase10 <- df_phase10$ph
df_kastner_model$phase10_CW <- df_phase10_CW$ph
df_kastner_model$phase10_CCW <- df_phase10_CCW$ph

cutData <- function( dataIn, subjId, roiInd, varThr ) {
  dataSel <- dataIn[ dataIn$participant==subjId & 
                       dataIn$ROI_name_basic==roiInd, ]
  varCut <- dataSel$var > quantile( dataSel$adjR2_05_CW_CCW_avg, varThr )
  dataSel <- dataSel[ varCut, ]
  return( dataSel )
}

##### for each participant and roi, cut the data #####
flagDf <- 1
for ( participantId in unique( df_kastner_model$participant ) ) { # participantId <- 'ASM16'
  df_filt_subj <- df_kastner_model[ df_kastner_model$participant==participantId, ]
  for ( roiId in unique( df_kastner_model$ROI_name_basic ) ) { # roiId <- 'V1'
    dataDfTemp <- cutData( df_kastner_model, participantId, roiId, 0.6  ) 
    if (flagDf==1) { df_kastner_model_filt_clean <- dataDfTemp; flagDf <- 0  }
    if (flagDf==0) { df_kastner_model_filt_clean <- rbind( df_kastner_model_filt_clean, dataDfTemp ) }
  }
}
str( df_kastner_model_filt_clean ) 
table( df_kastner_model_filt_clean$ROI_name_basic, df_kastner_model_filt_clean$participant )
df_kastner_model_filt_clean$ROI_name_basic <- as.character( df_kastner_model_filt_clean$ROI_name_basic )
df_kastner_model_filt_clean <- subset( df_kastner_model_filt_clean, ROI_name_basic=='V1' | ROI_name_basic=='V2' | ROI_name_basic=='V3'  | ROI_name_basic=='hV4' 
                                       | ROI_name_basic=='MST' | ROI_name_basic=='hMT' | ROI_name_basic=='V3a' | ROI_name_basic=='V3b' | ROI_name_basic=='IPS0' | ROI_name_basic=='IPS1'
                                       | ROI_name_basic=='IPS2' | ROI_name_basic=='IPS3' | ROI_name_basic=='IPS4' )
df_kastner_model_filt_clean$ROI_name_basic <- ifelse( df_kastner_model_filt_clean$ROI_name_basic == 'hMT' | df_kastner_model_filt_clean$ROI_name_basic == 'MST', 'hMT/MST', df_kastner_model_filt_clean$ROI_name_basic )
roiLevels <- unique(  df_kastner_model_filt_clean$ROI_name_basic )
df_kastner_model_filt_clean$roi_sorted <- factor( df_kastner_model_filt_clean$ROI_name_basic, levels=roiLevels ) 

##### all rois, kastner var explained and analysis, model 05, average #####
x11(width=4.5, height=4.5)
str( df_kastner_model )
df_kastner_model_aggregate_var <- aggregate( adjR2_05 ~ participant * roi_sorted, mean, data= df_kastner_model_filt_clean )
par(mfrow=c(1,1))
bp <- boxplot( df_kastner_model_aggregate_var$adjR2_05 ~ df_kastner_model_aggregate_var$roi_sorted, ylim=c(0,1), frame=FALSE, las=2,
         col=subjCols, outline=FALSE, xlab='ROI', ylab='adj. R^2 (proportion)', cex.lab=1.35, cex.axis=1.25 )
sMod <- summary( lmer( adjR2_05 ~ roi_sorted + ( 1 | participant ), data=df_kastner_model_aggregate_var ) )
for (i in 1:length(sMod$coefficients[,5])) { 
  if ( sMod$coefficients[i,5] < 0.001 ) {
    if (i==1) { text( i, bp$stats[5,i]+0.1, '*', cex=3, col='red' ) }
    if (i!=1) { text( i, bp$stats[5,i]+0.1, '*', cex=3 ) }
  }
}
dev.copy2pdf( file=sprintf('%s/%s/kastner05_varExp_adjR2_average.pdf', mainDir, figuresOutput ), width=4.5, height=4.5 )
graphics.off()

##### all rois, kastner var explained and analysis, model 05, average, quantiles #####

subsetDF_quantile <- c()
for ( p in unique( df_kastner_model_filt_clean$participant ) ) { # p <- 'ASM16'
  subsetDF <- subset( df_kastner_model_filt_clean, participant==p )
  subsetDF$adjR2_05_quantile <- as.numeric( cut( subsetDF$adjR2_05, 
                                               breaks=quantile( subsetDF$adjR2_05, seq( 0, 1, 0.01 ) ), 
                                               include.lowest = TRUE ) ) / 100
  subsetDF_quantile <- rbind( subsetDF_quantile, subsetDF )
}

x11(width=4.5, height=4.5)
str( subsetDF_quantile )
df_kastner_model_aggregate_var <- aggregate( adjR2_05_quantile ~ participant * roi_sorted, mean, data=subsetDF_quantile )
par(mfrow=c(1,1))
bp <- boxplot( df_kastner_model_aggregate_var$adjR2_05_quantile ~ df_kastner_model_aggregate_var$roi_sorted, ylim=c(0,1), frame=FALSE, las=2,
               col=subjCols, outline=FALSE, xlab='ROI', ylab='adj. R^2 (proportion)', cex.lab=1.35, cex.axis=1.25 )
sMod <- summary( lmer( adjR2_05_quantile ~ roi_sorted + ( 1 | participant ), data=df_kastner_model_aggregate_var ) )
for (i in 1:length(sMod$coefficients[,5])) { 
  if ( sMod$coefficients[i,5] < 0.001 ) {
    if (i==1) { text( i, bp$stats[5,i]+0.1, '*', cex=3, col='red' ) }
    if (i!=1) { text( i, bp$stats[5,i]+0.1, '*', cex=3 ) }
  }
}
dev.copy2pdf( file=sprintf('%s/%s/kastner05_varExp_adjR2_average_quantiles.pdf', mainDir, figuresOutput ), width=4.5, height=4.5 )
graphics.off()

##### all rois, kastner var explained and analysis, model 05, CW #####
x11(width=4.5, height=4.5)
str( df_kastner_model )
df_kastner_model_aggregate_var <- aggregate( adjR2_05_CW ~ participant * roi_sorted, mean, data= df_kastner_model_filt_clean )
par(mfrow=c(1,1))
bp <- boxplot( df_kastner_model_aggregate_var$adjR2_05_CW ~ df_kastner_model_aggregate_var$roi_sorted, ylim=c(0,1), frame=FALSE, las=2,
               col=subjCols, outline=FALSE, xlab='ROI', ylab='adj. R^2 (proportion)', cex.lab=1.35, cex.axis=1.25 )
sMod <- summary( lmer( adjR2_05_CW ~ roi_sorted + ( 1 | participant ), data=df_kastner_model_aggregate_var ) )
for (i in 1:length(sMod$coefficients[,5])) { 
  if ( sMod$coefficients[i,5] < 0.001 ) {
    if (i==1) { text( i, bp$stats[5,i]+0.1, '*', cex=3, col='red' ) }
    if (i!=1) { text( i, bp$stats[5,i]+0.1, '*', cex=3 ) }
  }
}
dev.copy2pdf( file=sprintf('%s/%s/kastner05_varExp_adjR2_CW.pdf', mainDir, figuresOutput ), width=4.5, height=4.5 )
graphics.off()

##### all rois, kastner var explained and analysis, model 05, CCW #####
x11(width=4.5, height=4.5)
str( df_kastner_model )
df_kastner_model_aggregate_var <- aggregate( adjR2_05_CCW ~ participant * roi_sorted, mean, data= df_kastner_model_filt_clean )
par(mfrow=c(1,1))
bp <- boxplot( df_kastner_model_aggregate_var$adjR2_05_CCW ~ df_kastner_model_aggregate_var$roi_sorted, ylim=c(0,1), frame=FALSE, las=2,
               col=subjCols, outline=FALSE, xlab='ROI', ylab='adj. R^2 (proportion)', cex.lab=1.35, cex.axis=1.25 )
sMod <- summary( lmer( adjR2_05_CCW ~ roi_sorted + ( 1 | participant ), data=df_kastner_model_aggregate_var ) )
for (i in 1:length(sMod$coefficients[,5])) { 
  if ( sMod$coefficients[i,5] < 0.001 ) {
    if (i==1) { text( i, bp$stats[5,i]+0.1, '*', cex=3, col='red' ) }
    if (i!=1) { text( i, bp$stats[5,i]+0.1, '*', cex=3 ) }
  }
}
dev.copy2pdf( file=sprintf('%s/%s/kastner05_varExp_adjR2_CCW.pdf', mainDir, figuresOutput ), width=4.5, height=4.5 )
graphics.off()

##### all rois, kappa par, model 05, average #####
x11(width=4.5, height=4.5)
df_kastner_model_aggregate_kappa <- aggregate( kappa ~ participant * roi_sorted, mean, data= df_kastner_model_filt_clean ) 
par(mfrow=c(1,1))
bp <- boxplot( df_kastner_model_aggregate_kappa$kappa ~ df_kastner_model_aggregate_kappa$roi_sorted, ylim=c(0.2,1), frame=FALSE, las=2,
         col=subjCols, outline=FALSE, xlab='ROI', ylab='kappa par. (tuning width)', cex.lab=1.35, cex.axis=1.25 )
sMod <- summary( lmer( kappa ~ roi_sorted + ( 1 | participant ), data=df_kastner_model_aggregate_kappa ) )
for (i in 1:length(sMod$coefficients[,5])) { 
  if ( sMod$coefficients[i,5] < 0.001 ) {
    if (i==1) { text( i, bp$stats[5,i]+0.1, '*', cex=3, col='red' ) }
    if (i!=1) { text( i, bp$stats[5,i]+0.1, '*', cex=3 ) }
  }
}
dev.copy2pdf( file=sprintf('%s/%s/kastner05_kappaPar_average.pdf', mainDir, figuresOutput ), width=4.5, height=4.5 )
graphics.off()

##### all rois, kappa par, model 05, CW #####
x11(width=4.5, height=4.5)
df_kastner_model_aggregate_kappa <- aggregate( kappa05_CW ~ participant * roi_sorted, mean, data= df_kastner_model_filt_clean ) 
par(mfrow=c(1,1))
bp <- boxplot( df_kastner_model_aggregate_kappa$kappa05_CW ~ df_kastner_model_aggregate_kappa$roi_sorted, ylim=c(0.2,1), frame=FALSE, las=2,
               col=subjCols, outline=FALSE, xlab='ROI', ylab='kappa par. (tuning width)', cex.lab=1.35, cex.axis=1.25 )
sMod <- summary( lmer( kappa05_CW ~ roi_sorted + ( 1 | participant ), data=df_kastner_model_aggregate_kappa ) )
for (i in 1:length(sMod$coefficients[,5])) { 
  if ( sMod$coefficients[i,5] < 0.001 ) {
    if (i==1) { text( i, bp$stats[5,i]+0.1, '*', cex=3, col='red' ) }
    if (i!=1) { text( i, bp$stats[5,i]+0.1, '*', cex=3 ) }
  }
}
dev.copy2pdf( file=sprintf('%s/%s/kastner05_kappaPar_CW.pdf', mainDir, figuresOutput ), width=4.5, height=4.5 )
graphics.off()

##### all rois, kappa par, model 05, CCW #####
x11(width=4.5, height=4.5)
df_kastner_model_aggregate_kappa <- aggregate( kappa05_CCW ~ participant * roi_sorted, mean, data= df_kastner_model_filt_clean ) 
par(mfrow=c(1,1))
bp <- boxplot( df_kastner_model_aggregate_kappa$kappa05_CCW ~ df_kastner_model_aggregate_kappa$roi_sorted, ylim=c(0.2,1), frame=FALSE, las=2,
               col=subjCols, outline=FALSE, xlab='ROI', ylab='kappa par. (tuning width)', cex.lab=1.35, cex.axis=1.25 )
sMod <- summary( lmer( kappa05_CCW ~ roi_sorted + ( 1 | participant ), data=df_kastner_model_aggregate_kappa ) )
for (i in 1:length(sMod$coefficients[,5])) { 
  if ( sMod$coefficients[i,5] < 0.001 ) {
    if (i==1) { text( i, bp$stats[5,i]+0.1, '*', cex=3, col='red' ) }
    if (i!=1) { text( i, bp$stats[5,i]+0.1, '*', cex=3 ) }
  }
}
dev.copy2pdf( file=sprintf('%s/%s/kastner05_kappaPar_CCW.pdf', mainDir, figuresOutput ), width=4.5, height=4.5 )
graphics.off()

##### all rois, model comparison #####
df_kastner_model_aggregate_adjR2 <- aggregate( adjR2_05 ~ participant * roi_sorted, mean, data= df_kastner_model_filt_clean ) 
df_kastner_model_aggregate_adjR2_10 <- aggregate( adjR2_10 ~ participant * roi_sorted, mean, data= df_kastner_model_filt_clean ) 
df_kastner_model_aggregate_adjR2_both <- aggregate( adjR2_both ~ participant * roi_sorted, mean, data= df_kastner_model_filt_clean ) 
df_kastner_model_aggregate_adjR2_phase05 <- aggregate( adjR2_phase05 ~ participant * roi_sorted, mean, data= df_kastner_model_filt_clean ) 
df_kastner_model_aggregate_adjR2_phase10 <- aggregate( adjR2_phase10 ~ participant * roi_sorted, mean, data= df_kastner_model_filt_clean ) 

df_kastner_model_aggregate_adjR2$adjR2_10 <- df_kastner_model_aggregate_adjR2_10$adjR2_10
df_kastner_model_aggregate_adjR2$adjR2_both <- df_kastner_model_aggregate_adjR2_both$adjR2_both
df_kastner_model_aggregate_adjR2$adjR2_phase05 <- df_kastner_model_aggregate_adjR2_phase05$adjR2_phase05
df_kastner_model_aggregate_adjR2$adjR2_phase10 <- df_kastner_model_aggregate_adjR2_phase10$adjR2_phase10

x11(width=4.5, height=4.5)
par(mfrow=c(1,1))
boxplot( df_kastner_model_aggregate_adjR2$adjR2_05 ~ df_kastner_model_aggregate_adjR2$roi_sorted, 
         add=FALSE, boxwex=0.3, ylim=c(0, 0.6), at = 1:12, xlim=c(0,13), yaxs='i', axes=TRUE, outline=FALSE, lty=1, col='gray50',
         frame=FALSE, las=2, cex.lab=1.35, cex.axis=1.25, xlab='ROI', ylab='adj. R^2 (proportion)' )
bp <- boxplot( df_kastner_model_aggregate_adjR2$adjR2_phase05 ~ df_kastner_model_aggregate_adjR2$roi_sorted, 
         add=TRUE, boxwex=0.3, at = 1:12+0.3, xlim=c(0,13), yaxs='i', frame=FALSE, axes=FALSE, outline=FALSE, lty=1, col='red' )
for (i in 1:dim(bp$stats)[2]) {
  selectedLevel <- levels( df_kastner_model_aggregate_adjR2$roi_sorted )[i]
  t.out <- t.test( df_kastner_model_aggregate_adjR2$adjR2_phase05[ df_kastner_model_aggregate_adjR2$roi_sorted==selectedLevel ],
          df_kastner_model_aggregate_adjR2$adjR2_05[ df_kastner_model_aggregate_adjR2$roi_sorted==selectedLevel ], paire=TRUE )
  if ( t.out$p.value < 0.01/14 ) {
    text( i, bp$stats[5,i]+0.2, '*', cex=3 )
  }
}
legend( 6.8, 0.12, c('model','phase encoded'), bty='n', fill=c('gray50','red') )
dev.copy2pdf( file=sprintf('%s/%s/kastner05_modelComparison.pdf', mainDir, figuresOutput ), width=4.5, height=4.5 )
graphics.off()


##### full comparison (phase_05, phase_10, model05, model10 etc etc), separate manuscript #####

x11(width=9, height=4.5)
par(mfrow=c(1,1))
boxplot( df_kastner_model_aggregate_adjR2$adjR2_both ~ df_kastner_model_aggregate_adjR2$roi_sorted,
         add=FALSE, boxwex=0.2, ylim=c(0, 0.6), at = 1:12, xlim=c(0,13), yaxs='i', axes=TRUE, outline=FALSE, lty=1, col='gray50',
         frame=FALSE, las=2, cex.lab=1.35, cex.axis=1.25, xlab='ROI', ylab='adj. R^2 (proportion)' )
boxplot( df_kastner_model_aggregate_adjR2$adjR2_05 ~ df_kastner_model_aggregate_adjR2$roi_sorted,
         add=TRUE, frame=FALSE, las=1, boxwex=0.2, at = 1:12-0.25, xlim=c(0,13), yaxs='i', outline=FALSE, lty=1, col='lightblue', axes=FALSE )
boxplot( df_kastner_model_aggregate_adjR2$adjR2_phase05 ~ df_kastner_model_aggregate_adjR2$roi_sorted,
         add=TRUE, boxwex=0.2, at = 1:12+0.25, xlim=c(0,13), yaxs='i', frame=FALSE, axes=FALSE, outline=FALSE, lty=1, col='red' )
boxplot( df_kastner_model_aggregate_adjR2$adjR2_10 ~ df_kastner_model_aggregate_adjR2$roi_sorted,
         add=TRUE, boxwex=0.2, at = 1:12-0.06, xlim=c(0,14), yaxs='i', frame=FALSE, axes=FALSE, outline=FALSE, lty=1, col='orange' )
boxplot( df_kastner_model_aggregate_adjR2$adjR2_phase10 ~ df_kastner_model_aggregate_adjR2$roi_sorted,
         add=TRUE, boxwex=0.2, at = 1:12-0.25-0.06, xlim=c(0,13), yaxs='i', frame=FALSE, axes=FALSE, outline=FALSE, lty=1, col='blue' )
abline( v=seq(1.5,11.5), lwd=0.5, lty=2 )
points( c(0,2), c(0.55, 0.55),  pch=15, cex=1.5, col=c('gray50', 'lightblue') )
points( c(0,2), c(0.50, 0.50),  pch=15, cex=1.5, col=c('red','orange') )
points( c(0), c(0.45),  pch=15, cex=1.5, col=c('blue') )
text( 0.5, 0.55, 'both'); text( 2.9, 0.55, 'model05'); text( 0.8, 0.50, 'phase05'); text( 2.9, 0.50, 'model10'); text( 0.8, 0.45, 'phase10');
dev.copy2pdf( file=sprintf('%s/%s/kastner_full_modelComparison.pdf', mainDir, figuresOutput ), width=9, height=4.5 )
graphics.off()

##### all rois, both multiplicative par, separate manuscript #####
x11(width=4.5, height=4.5)
df_kastner_model_aggregate_both_ratio <- aggregate( both_ratio ~ participant * roi_sorted, mean, data=df_kastner_model_filt_clean ) 
par(mfrow=c(1,1))
boxplot( df_kastner_model_aggregate_both_ratio$both_ratio ~ df_kastner_model_aggregate_both_ratio$roi_sorted, ylim=c(0,1),
         frame=FALSE, col=subjCols, las=2, cex.lab=1.35, cex.axis=1.25, xlab='ROI', ylab='05/10 ratio', outline=FALSE )
summary( lm( both_ratio ~ roi_sorted, data=df_kastner_model_aggregate_both_ratio) )
dev.copy2pdf( file=sprintf('%s/%s/kastner05_10_ratio.pdf', mainDir, figuresOutput ), width=4.5, height=4.5 )
graphics.off()

##### all rois, both multiplicative par, separate manuscript, rescaled #####
x11(width=4.5, height=4.5)
df_kastner_model_aggregate_both_ratio <- aggregate( both_ratio ~ participant * roi_sorted, mean, data=df_kastner_model_filt_clean ) 
par(mfrow=c(1,1))
boxplot( df_kastner_model_aggregate_both_ratio$both_ratio ~ df_kastner_model_aggregate_both_ratio$roi_sorted, ylim=c(0.6,1),
         frame=FALSE, col=subjCols, las=2, cex.lab=1.35, cex.axis=1.25, xlab='ROI', ylab='05/10 ratio', outline=FALSE )
summary( lm( both_ratio ~ roi_sorted, data=df_kastner_model_aggregate_both_ratio) )
dev.copy2pdf( file=sprintf('%s/%s/kastner05_10_ratio_rescaled.pdf', mainDir, figuresOutput ), width=4.5, height=4.5 )
graphics.off()


##### with histograms, density, modelling05 #####
df_kastner_model_filt_clean$subject <- df_kastner_model_filt_clean$participant
df_kastner_model_filt_clean$roi_coarse <- ifelse( ( df_kastner_model_filt_clean$roi_sorted=='IPS0' | df_kastner_model_filt_clean$roi_sorted=='IPS1' ), 'IPS0/IPS1',
                                          ifelse( ( df_kastner_model_filt_clean$roi_sorted=='IPS2' | df_kastner_model_filt_clean$roi_sorted=='IPS3' ), 'IPS2/IPS3',
                                          ifelse( ( df_kastner_model_filt_clean$roi_sorted=='IPS4' ), 'IPS4',
                                          ifelse( ( df_kastner_model_filt_clean$roi_sorted=='MST' | df_kastner_model_filt_clean$roi_sorted=='hMT' ), 'MST/hMT',
                                          ifelse( ( df_kastner_model_filt_clean$roi_sorted=='V3a' | df_kastner_model_filt_clean$roi_sorted=='V3b' ), 'V3a/V3b',        
                                          as.character( df_kastner_model_filt_clean$roi_sorted ) ) ) ) ) )
df_kastner_model_filt_clean$roi_coarse <- factor( df_kastner_model_filt_clean$roi_coarse, levels=unique( df_kastner_model_filt_clean$roi_coarse ) )

df_kastner_model_filt_clean$side <- ifelse( df_kastner_model_filt_clean$hemiName=='lh', 'left', 'right' )


#### ROC on average data model05, plot function ####

plotDirectionDistribution <- function( df_kastner_model_filt_clean, outputSuffix ) {

  nBreaks <- 10
  dfFlag <- 1
  
  for ( subjSel in unique( df_kastner_model_filt_clean$subject )  ) { # subjSel <- 'ASM16'
    x11( width=15, height=2 )
    par( mfrow=c(1,9) )
    for ( roiSel in levels( df_kastner_model_filt_clean$roi_coarse )  ) { # roiSel <- 'V1'
      
      df_kastner_model_filt_subj_left <- subset( df_kastner_model_filt_clean, subject==subjSel & roi_coarse==roiSel & side=='left' )
      ds_left <- hist( df_kastner_model_filt_subj_left$dv, plot=FALSE, breaks = seq(0,2*pi, length.out=nBreaks) )
      df_kastner_model_filt_subj_right <- subset( df_kastner_model_filt_clean, subject==subjSel & roi_coarse==roiSel & side=='right' )
      ds_right <- hist( df_kastner_model_filt_subj_right$dv, plot=FALSE, breaks = seq(0,2*pi, length.out=nBreaks) )
      
      # density
      maxD <- max( c( max(ds_left$density), max(ds_right$density) )  )
      plot( 0, type='n', ylab='density', xlab='mu (rad)', main=roiSel, xlim=c(0,2*pi), ylim=c(0,maxD), axes=FALSE, bty='n' )
      lines( ds_left$mids, ds_left$density , col='darkorange', lwd=3, lty=1 )
      lines( ds_right$mids, ds_right$density, col='darkblue', lwd=3, lty=1 )
      axis( 1, c(0,pi,2*pi), round( c(0,pi,2*pi),2 ) ); axis(2, c(0,round(maxD,2)), las=1 )
      
      subjArrayLeft <-  rep( subjSel, length( ds_left$density ) )
      sideArrayLeft <- rep( 'left', length( ds_left$density ) )
      roiArrayLeft <-  rep( roiSel, length( ds_left$density ) )
      densityArrayLeft <-  ds_left$density
      midsArrayLeft <- ds_left$mids
      
      subjArrayRight <-  rep( subjSel, length( ds_right$density ) )
      sideArrayRight <- rep( 'right', length( ds_right$density ) )
      roiArrayRight <-  rep( roiSel, length( ds_right$density ) )
      densityArrayRight <-  ds_right$density
      midsArrayRight <- ds_right$mids
      
      subjArray <- c( subjArrayLeft, subjArrayRight )
      sideArray <- c( sideArrayLeft, sideArrayRight )
      roiArray <- c( roiArrayLeft, roiArrayRight )
      densityArray <- c( densityArrayLeft, densityArrayRight )
      midsArray <- c( midsArrayLeft, midsArrayRight )
      
      subjDF <- data.frame( subjArray, sideArray, roiArray, densityArray, midsArray  )
      
      if (dfFlag == 1) { subjDFOutDensity <- subjDF; dfFlag <- 2 }
      if (dfFlag>1) { subjDFOutDensity <- rbind( subjDFOutDensity, subjDF ) } 
    }
    dev.copy2pdf( file=sprintf('%s/%s/directionResultsModelling05/subject_%s_%s.pdf', mainDir, figuresOutput, subjSel, outputSuffix ), width=15, height=2 )
    graphics.off()
  }
  
  aggDFDensity <- aggregate( densityArray ~ midsArray*sideArray*roiArray, mean, data=subjDFOutDensity )
  aggDFDensitySd <- aggregate( densityArray ~ midsArray*sideArray*roiArray, sd, data=subjDFOutDensity )
  aggDFDensitySd$Se <- aggDFDensitySd$densityArray  / sqrt( length( unique( df_kastner_model_filt_clean$subject ) ) )
  aggDFDensity$roiArray <- as.factor( aggDFDensity$roiArray )
  levels( aggDFDensity$roiArray ) <- levels( aggDFDensity$roiArray )[c(6,7,8,9,5,1,2,3,4)]
  
  x11( width=15, height=2 )
  par( mfrow=c(1,9) )
  for ( roiIdx in levels( aggDFDensity$roiArray ) ) {
    selectedDf <- subset( aggDFDensity, roiArray==roiIdx )
    selectedDfSd <- subset( aggDFDensitySd, roiArray==roiIdx )
    maxD <- max( c( max(selectedDf$densityArray), max(selectedDf$densityArray) )  )
    maxD <- maxD + maxD*0.15
    plot( 0, type='n', ylab='density', xlab='mu (rad)', main=roiIdx, xlim=c(0,2*pi), ylim=c(0,maxD), axes=FALSE, bty='n' )
    lines( selectedDf$midsArray[ selectedDf$sideArray=='left' ], selectedDf$densityArray[ selectedDf$sideArray=='left' ], col='darkorange', lwd=3, lty=1 )
    lines( selectedDf$midsArray[ selectedDf$sideArray=='left' ], selectedDf$densityArray[ selectedDf$sideArray=='left' ]+selectedDfSd$Se[ selectedDf$sideArray=='left' ], col='orange', lwd=1, lty=1 )
    lines( selectedDf$midsArray[ selectedDf$sideArray=='left' ], selectedDf$densityArray[ selectedDf$sideArray=='left' ]-selectedDfSd$Se[ selectedDf$sideArray=='left' ], col='orange', lwd=1, lty=1 )
    
    lines( selectedDf$midsArray[ selectedDf$sideArray=='right' ], selectedDf$densityArray[ selectedDf$sideArray=='right' ], col='darkblue', lwd=3, lty=1 )
    lines( selectedDf$midsArray[ selectedDf$sideArray=='right' ], selectedDf$densityArray[ selectedDf$sideArray=='right' ]+selectedDfSd$Se[ selectedDf$sideArray=='right' ], col='darkblue', lwd=1, lty=1 )
    lines( selectedDf$midsArray[ selectedDf$sideArray=='right' ], selectedDf$densityArray[ selectedDf$sideArray=='right' ]-selectedDfSd$Se[ selectedDf$sideArray=='right' ], col='darkblue', lwd=1, lty=1 )
    axis( 1, c(0,pi,2*pi), round( c(0,pi,2*pi),2 ) ); axis(2, c(0,round(maxD,2)), las=1 )
  }
  dev.copy2pdf( file=sprintf('%s/%s/directionResultsModelling05/allSubjects_%s.pdf', mainDir, figuresOutput, outputSuffix ), width=15, height=2 )
  
  graphics.off()
  
  return( subjDFOutDensity )
  
}

#### ROC on average data model05 ####

df_kastner_model_filt_clean$dv <- atan2(
  ( sin( df_kastner_model_filt_clean$mu05_CW-pi ) + sin( df_kastner_model_filt_clean$mu05_CCW-pi ) ) / 2, 
  ( cos( df_kastner_model_filt_clean$mu05_CW-pi ) + cos( df_kastner_model_filt_clean$mu05_CCW-pi ) ) / 2
) + pi
subjDf_average <- plotDirectionDistribution( df_kastner_model_filt_clean, 'average' )

#### example plot
pROC_obj <- roc(df_kastner_model_filt_clean$side, df_kastner_model_filt_clean$dv,
                smoothed = TRUE,
                ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                print.auc=TRUE, show.thres=TRUE)
dev.copy2pdf( file=sprintf('%s/%s/exampleROC_curve.pdf', mainDir, figuresOutput ), width=4.5, height=4.5 )
graphics.off()

sOut <- c()
rOut <- c()
aucOut_average <- c()
df_kastner_model_filt_clean$participant <- as.factor( df_kastner_model_filt_clean$participant )
for ( s in levels( df_kastner_model_filt_clean$participant ) ) {
  for ( r in levels( df_kastner_model_filt_clean$roi_coarse )[1:9] ) {
    print( sprintf( 'subject: %s, roi: %s', s, r ) )
    dfLoop <- subset( df_kastner_model_filt_clean, participant==s & roi_coarse==r )
    aTemp <- roc( dfLoop$side, dfLoop$dv, print.auc=FALSE, print.auc = FALSE )
    auc <- round( as.numeric( aTemp$auc ), 3 )
    sOut <- c( sOut, s )
    rOut <- c( rOut, r )
    aucOut_average <- c( aucOut_average, auc )
  }
}

#rOut <- factor( rOut, levels=unique(rOut) )
#sOut <- factor( sOut, levels=unique(sOut) )
#x11( width=6, height=4.5 )
#boxplot( aucOut ~ rOut, frame=FALSE, las=2, ylab='AUC', xlab='ROI', ylim=c(0.2,1),
#         outline=FALSE, cex.lab=1.5, cex.axis=1.5, lwd=2 )
#summary( lmer( aucOut ~ rOut + ( 1 | sOut ) ) )
#dev.copy2pdf( file=sprintf('%s/%s/aucBoxPlot_curve_average.pdf', mainDir, figuresOutput ), width=6, height=4.5 )
#graphics.off()


#### ROC on CCW data model05 ####

df_kastner_model_filt_clean$dv <- df_kastner_model_filt_clean$mu05_CCW
plotDirectionDistribution( df_kastner_model_filt_clean, 'CCW' )

sOut <- c()
rOut <- c()
aucOut_CCW <- c()
df_kastner_model_filt_clean$participant <- as.factor( df_kastner_model_filt_clean$participant )
for ( s in levels( df_kastner_model_filt_clean$participant ) ) {
  for ( r in levels( df_kastner_model_filt_clean$roi_coarse )[1:9] ) {
    print( sprintf( 'subject: %s, roi: %s', s, r ) )
    dfLoop <- subset( df_kastner_model_filt_clean, participant==s & roi_coarse==r )
    aTemp <- roc( dfLoop$side, dfLoop$dv, print.auc=FALSE, print.auc = FALSE )
    auc <- round( as.numeric( aTemp$auc ), 3 )
    sOut <- c( sOut, s )
    rOut <- c( rOut, r )
    aucOut_CCW <- c( aucOut_CCW, auc )
  }
}

#### ROC on CW data model 05

df_kastner_model_filt_clean$dv <- df_kastner_model_filt_clean$mu05_CW
plotDirectionDistribution( df_kastner_model_filt_clean, 'CW' )

sOut <- c()
rOut <- c()
aucOut_CW <- c()
df_kastner_model_filt_clean$participant <- as.factor( df_kastner_model_filt_clean$participant )
for ( s in levels( df_kastner_model_filt_clean$participant ) ) {
  for ( r in levels( df_kastner_model_filt_clean$roi_coarse )[1:9] ) {
    print( sprintf( 'subject: %s, roi: %s', s, r ) )
    dfLoop <- subset( df_kastner_model_filt_clean, participant==s & roi_coarse==r )
    aTemp <- roc( dfLoop$side, dfLoop$dv, print.auc=FALSE, print.auc = FALSE )
    auc <- round( as.numeric( aTemp$auc ), 3 )
    sOut <- c( sOut, s )
    rOut <- c( rOut, r )
    aucOut_CW <- c( aucOut_CW, auc )
  }
}

rOut <- factor( rOut, levels=unique(rOut) )
sOut <- factor( sOut, levels=unique(sOut) )
x11( width=6, height=4.5 )
boxplot( aucOut_average ~ rOut, frame=FALSE, las=2, ylab='AUC', xlab='ROI', ylim=c(0.2,1),
         outline=FALSE, cex.lab=1.5, cex.axis=1.5, lwd=2, boxwex=0.2, at=1:9-0.4, col='gray25' )
boxplot( aucOut_CW ~ rOut,  frame=FALSE, axes=FALSE, outline=FALSE, las=2, ylab='', xlab='',
         lwd=2, boxwex=0.2, at=1:9-0.2, col='gray70', add=TRUE )
boxplot( aucOut_CCW ~ rOut,  frame=FALSE, axes=FALSE, outline=FALSE, las=2, ylab='', xlab='',
         lwd=2, boxwex=0.2, at=1:9, col='gray90', add=TRUE )
legend( 5, 0.6, col=c('gray25','gray70','gray90'), lwd=7, legend=c('average','CW','CCW'), cex=1.5, bty='n' )
summary( lmer( aucOut_average ~ rOut + ( 1 | sOut ) ) )
summary( lmer( aucOut_CW ~ rOut + ( 1 | sOut ) ) )
summary( lmer( aucOut_CCW ~ rOut + ( 1 | sOut ) ) )
dev.copy2pdf( file=sprintf('%s/%s/aucBoxPlot_curve_model05.pdf', mainDir, figuresOutput ), width=6, height=4.5 )
graphics.off()


##### with histograms, density, phase05 #####

plotDirectionDistribution <- function( df_kastner_model_filt_clean, outputSuffix ) {
  
  nBreaks <- 10
  dfFlag <- 1
  
  for ( subjSel in unique( df_kastner_model_filt_clean$participant )  ) {
    x11( width=15, height=2 )
    par( mfrow=c(1,9) )
    for ( roiSel in levels( df_kastner_model_filt_clean$roi_coarse )  ) {
      
      df_kastner_model_filt_subj_left <- subset( df_kastner_model_filt_clean, participant==subjSel & roi_coarse==roiSel & side=='left' )
      ds_left <- hist( df_kastner_model_filt_subj_left$dv, plot=FALSE, breaks = seq(0,2*pi, length.out=nBreaks) )
      df_kastner_model_filt_subj_right <- subset( df_kastner_model_filt_clean, participant==subjSel & roi_coarse==roiSel & side=='right' )
      ds_right <- hist( df_kastner_model_filt_subj_right$dv, plot=FALSE, breaks = seq(0,2*pi, length.out=nBreaks) )

      # density
      maxD <- max( c( max(ds_left$density), max(ds_right$density) )  )
      plot( 0, type='n', ylab='density', xlab='mu (rad)', main=roiSel, xlim=c(0,2*pi), ylim=c(0,maxD), axes=FALSE, bty='n' )
      lines( ds_left$mids, ds_left$density , col='darkorange', lwd=3, lty=1 )
      lines( ds_right$mids, ds_right$density, col='darkblue', lwd=3, lty=1 )
      axis( 1, c(0,pi,2*pi), round( c(0,pi,2*pi),2 ) ); axis(2, c(0,round(maxD,2)), las=1 )
      
      subjArrayLeft <-  rep( subjSel, length( ds_left$density ) )
      sideArrayLeft <- rep( 'left', length( ds_left$density ) )
      roiArrayLeft <-  rep( roiSel, length( ds_left$density ) )
      densityArrayLeft <-  ds_left$density
      midsArrayLeft <- ds_left$mids
      
      subjArrayRight <-  rep( subjSel, length( ds_right$density ) )
      sideArrayRight <- rep( 'right', length( ds_right$density ) )
      roiArrayRight <-  rep( roiSel, length( ds_right$density ) )
      densityArrayRight <-  ds_right$density
      midsArrayRight <- ds_right$mids
      
      subjArray <- c( subjArrayLeft, subjArrayRight )
      sideArray <- c( sideArrayLeft, sideArrayRight )
      roiArray <- c( roiArrayLeft, roiArrayRight )
      densityArray <- c( densityArrayLeft, densityArrayRight )
      midsArray <- c( midsArrayLeft, midsArrayRight )
      
      subjDF <- data.frame( subjArray, sideArray, roiArray, densityArray, midsArray  )
      
      if (dfFlag == 1) { subjDFOutDensity <- subjDF; dfFlag <- 2 }
      if (dfFlag>1) { subjDFOutDensity <- rbind( subjDFOutDensity, subjDF ) } 
    }
    dev.copy2pdf( file=sprintf('%s/%s/directionResultsPhase05/subject_%s_%s.pdf', mainDir, figuresOutput, subjSel, outputSuffix ), width=15, height=2 )
    graphics.off()
  }
  
  aggDFDensity <- aggregate( densityArray ~ midsArray*sideArray*roiArray, mean, data=subjDFOutDensity )
  aggDFDensitySd <- aggregate( densityArray ~ midsArray*sideArray*roiArray, sd, data=subjDFOutDensity )
  aggDFDensitySd$Se <- aggDFDensitySd$densityArray  / sqrt( length( unique( df_kastner_model_filt_clean$subject ) ) )
  aggDFDensity$roiArray <- as.factor( aggDFDensity$roiArray )
  levels( aggDFDensity$roiArray ) <- levels( aggDFDensity$roiArray )[c(6,7,8,9,5,1,2,3,4)]
  
  x11( width=15, height=2 )
  par( mfrow=c(1,9) )
  for ( roiIdx in levels( aggDFDensity$roiArray ) ) {
    selectedDf <- subset( aggDFDensity, roiArray==roiIdx )
    selectedDfSd <- subset( aggDFDensitySd, roiArray==roiIdx )
    maxD <- max( c( max(selectedDf$densityArray), max(selectedDf$densityArray) )  )
    maxD <- maxD + maxD*0.15
    
    plot( 0, type='n', ylab='density', xlab='mu (rad)', main=roiIdx, xlim=c(0,2*pi), ylim=c(0,maxD), axes=FALSE, bty='n' )
    lines( selectedDf$midsArray[ selectedDf$sideArray=='left' ], selectedDf$densityArray[ selectedDf$sideArray=='left' ], col='darkorange', lwd=3, lty=1 )
    lines( selectedDf$midsArray[ selectedDf$sideArray=='left' ], selectedDf$densityArray[ selectedDf$sideArray=='left' ]+selectedDfSd$Se[ selectedDf$sideArray=='left' ], col='orange', lwd=1, lty=1 )
    lines( selectedDf$midsArray[ selectedDf$sideArray=='left' ], selectedDf$densityArray[ selectedDf$sideArray=='left' ]-selectedDfSd$Se[ selectedDf$sideArray=='left' ], col='orange', lwd=1, lty=1 )
    
    lines( selectedDf$midsArray[ selectedDf$sideArray=='right' ], selectedDf$densityArray[ selectedDf$sideArray=='right' ], col='darkblue', lwd=3, lty=1 )
    lines( selectedDf$midsArray[ selectedDf$sideArray=='right' ], selectedDf$densityArray[ selectedDf$sideArray=='right' ]+selectedDfSd$Se[ selectedDf$sideArray=='right' ], col='darkblue', lwd=1, lty=1 )
    lines( selectedDf$midsArray[ selectedDf$sideArray=='right' ], selectedDf$densityArray[ selectedDf$sideArray=='right' ]-selectedDfSd$Se[ selectedDf$sideArray=='right' ], col='darkblue', lwd=1, lty=1 )
    axis( 1, c(0,pi,2*pi), round( c(0,pi,2*pi),2 ) ); axis(2, c(0,round(maxD,2)), las=1 )
    
  }
  dev.copy2pdf( file=sprintf('%s/%s/directionResultsPhase05/allSubjects_%s.pdf', mainDir, figuresOutput, outputSuffix ), width=15, height=2 )
  graphics.off()
}

#### ROC on average data phase 05

dvCW <- ( (df_kastner_model_filt_clean$phase05_CW+2) %% (2*pi) ) * -1 + 2*pi
dvCCW <- (df_kastner_model_filt_clean$phase05_CCW+2) %% (2*pi)
df_kastner_model_filt_clean$dv <- atan2(
  ( sin( dvCW-pi ) + sin( dvCCW-pi ) ) / 2, 
  ( cos( dvCW-pi ) + cos( dvCCW-pi ) ) / 2
) + pi
str( df_kastner_model_filt_clean )
plotDirectionDistribution( df_kastner_model_filt_clean, 'average' )

sOut <- c()
rOut <- c()
aucOut_average <- c()
df_kastner_model_filt_clean$participant <- as.factor( df_kastner_model_filt_clean$participant )
for ( s in levels( df_kastner_model_filt_clean$participant ) ) {
  for ( r in levels( df_kastner_model_filt_clean$roi_coarse )[1:9] ) {
    print( sprintf( 'subject: %s, roi: %s', s, r ) )
    dfLoop <- subset( df_kastner_model_filt_clean, participant==s & roi_coarse==r )
    aTemp <- roc( dfLoop$side, dfLoop$dv, print.auc=FALSE, print.auc = FALSE )
    auc <- round( as.numeric( aTemp$auc ), 3 )
    sOut <- c( sOut, s )
    rOut <- c( rOut, r )
    aucOut_average <- c( aucOut_average, auc )
  }
}

#### ROC on CCW data phase 05

df_kastner_model_filt_clean$dv <- (df_kastner_model_filt_clean$phase05_CCW+2) %% (2*pi)
plotDirectionDistribution( df_kastner_model_filt_clean, 'CCW' )

sOut <- c()
rOut <- c()
aucOut_CCW <- c()
df_kastner_model_filt_clean$participant <- as.factor( df_kastner_model_filt_clean$participant )
for ( s in levels( df_kastner_model_filt_clean$participant ) ) {
  for ( r in levels( df_kastner_model_filt_clean$roi_coarse )[1:9] ) {
    print( sprintf( 'subject: %s, roi: %s', s, r ) )
    dfLoop <- subset( df_kastner_model_filt_clean, participant==s & roi_coarse==r )
    aTemp <- roc( dfLoop$side, dfLoop$dv, print.auc=FALSE, print.auc = FALSE )
    auc <- round( as.numeric( aTemp$auc ), 3 )
    sOut <- c( sOut, s )
    rOut <- c( rOut, r )
    aucOut_CCW <- c( aucOut_CCW, auc )
  }
}

#### ROC on CW data phase 05

df_kastner_model_filt_clean$dv <- ( (df_kastner_model_filt_clean$phase05_CW+2) %% (2*pi) ) * -1 + 2*pi
plotDirectionDistribution( df_kastner_model_filt_clean, 'CW' )

sOut <- c()
rOut <- c()
aucOut_CW <- c()
df_kastner_model_filt_clean$participant <- as.factor( df_kastner_model_filt_clean$participant )
for ( s in levels( df_kastner_model_filt_clean$participant ) ) {
  for ( r in levels( df_kastner_model_filt_clean$roi_coarse )[1:9] ) {
    print( sprintf( 'subject: %s, roi: %s', s, r ) )
    dfLoop <- subset( df_kastner_model_filt_clean, participant==s & roi_coarse==r )
    aTemp <- roc( dfLoop$side, dfLoop$dv, print.auc=FALSE, print.auc = FALSE )
    auc <- round( as.numeric( aTemp$auc ), 3 )
    sOut <- c( sOut, s )
    rOut <- c( rOut, r )
    aucOut_CW <- c( aucOut_CW, auc )
  }
}

rOut <- factor( rOut, levels=unique(rOut) )
sOut <- factor( sOut, levels=unique(sOut) )
x11( width=6, height=4.5 )
boxplot( aucOut_average ~ rOut, frame=FALSE, las=2, ylab='AUC', xlab='ROI', ylim=c(0.2,1),
         outline=FALSE, cex.lab=1.5, cex.axis=1.5, lwd=2, boxwex=0.2, at=1:9-0.4, col='gray25' )
boxplot( aucOut_CW ~ rOut,  frame=FALSE, axes=FALSE, outline=FALSE, las=2, ylab='', xlab='',
         lwd=2, boxwex=0.2, at=1:9-0.2, col='gray70', add=TRUE )
boxplot( aucOut_CCW ~ rOut,  frame=FALSE, axes=FALSE, outline=FALSE, las=2, ylab='', xlab='',
         lwd=2, boxwex=0.2, at=1:9, col='gray90', add=TRUE )
legend( 5, 0.6, col=c('gray25','gray70','gray90'), lwd=7, legend=c('average','CW','CCW'), cex=1.5, bty='n' )
summary( lmer( aucOut_average ~ rOut + ( 1 | sOut ) ) )
summary( lmer( aucOut_CW ~ rOut + ( 1 | sOut ) ) )
summary( lmer( aucOut_CCW ~ rOut + ( 1 | sOut ) ) )
dev.copy2pdf( file=sprintf('%s/%s/aucBoxPlot_curve_phase.pdf', mainDir, figuresOutput ), width=6, height=4.5 )
graphics.off()


str( df_kastner_model_filt_clean )
table( df_kastner_model_filt_clean$participant, df_kastner_model_filt_clean$roi_coarse )
df_kastner_model_filt_clean$dv <- atan2(
  ( sin( df_kastner_model_filt_clean$mu05_CW ) + sin( df_kastner_model_filt_clean$mu05_CCW ) ) / 2, 
  ( cos( df_kastner_model_filt_clean$mu05_CW ) + cos( df_kastner_model_filt_clean$mu05_CCW ) ) / 2
) + pi / 2
df_kastner_model_filt_clean$dv <- df_kastner_model_filt_clean$mu
tempAa <- subset( df_kastner_model_filt_clean, participant=='HMR28' & roi_coarse=='V3a/V3b' & side=='left' )
coordsPlot <- pol2cart( cbind( tempAa$dv, tempAa$kappa ) )
xPlot <- coordsPlot[,1]
yPlot <- coordsPlot[,2]
plot( xPlot, yPlot, cex=0.5, xlim=c(-1,1), ylim=c(-1,1), type='p' )
selectedPoints <- tempAa$dv > ( pi/4 ) & tempAa$dv < ( pi/2 + pi/4 )
points( xPlot[selectedPoints], yPlot[selectedPoints], col='red' )

polar( tempAa$dv, tempAa$kappa, 'n' )


#####

#par( mfrow=c(1,5) )
#for ( roiIdx in levels( aggDFDensity$roiArray ) ) {
  # selectedDf <- subset( aggDFDensity, roiArray==roiIdx )
  # selectedDfSd <- subset( aggDFDensitySd, roiArray==roiIdx )
  # maxD <- max( c( max(selectedDf$densityArray), max(selectedDf$densityArray) )  )
  # maxD <- maxD + maxD*0.15
  # plot( 0, type='n', ylab='density', xlab='mu (rad)', main=roiIdx, xlim=c(0,2*pi), ylim=c(0,maxD), axes=FALSE, bty='n' )
  # lines( selectedDf$midsArray[ selectedDf$sideArray=='left' ], selectedDf$densityArray[ selectedDf$sideArray=='left' ], col='darkorange', lwd=3, lty=1 )
  # lines( selectedDf$midsArray[ selectedDf$sideArray=='left' ], selectedDf$densityArray[ selectedDf$sideArray=='left' ]+selectedDfSd$Se[ selectedDf$sideArray=='left' ], col='orange', lwd=1, lty=1 )
  # lines( selectedDf$midsArray[ selectedDf$sideArray=='left' ], selectedDf$densityArray[ selectedDf$sideArray=='left' ]-selectedDfSd$Se[ selectedDf$sideArray=='left' ], col='orange', lwd=1, lty=1 )
  # 
  # lines( selectedDf$midsArray[ selectedDf$sideArray=='right' ], selectedDf$densityArray[ selectedDf$sideArray=='right' ], col='darkblue', lwd=3, lty=1 )
  # lines( selectedDf$midsArray[ selectedDf$sideArray=='right' ], selectedDf$densityArray[ selectedDf$sideArray=='right' ]+selectedDfSd$Se[ selectedDf$sideArray=='right' ], col='darkblue', lwd=1, lty=1 )
  # lines( selectedDf$midsArray[ selectedDf$sideArray=='right' ], selectedDf$densityArray[ selectedDf$sideArray=='right' ]-selectedDfSd$Se[ selectedDf$sideArray=='right' ], col='darkblue', lwd=1, lty=1 )
  # axis( 1, c(0,pi,2*pi), round( c(0,pi,2*pi),2 ) ); axis(2, c(0,round(maxD,2)), las=1 )




#########################################
#########################################
## model kastner (filtered on prf data)##
#########################################
#########################################
str(df_model05)
df_kastner_model <- df_model05
df_kastner_model$logLik_10 <- df_model10$logLik
df_kastner_model$AIC_10 <- df_model10$AIC
df_kastner_model$adjR2_10 <- df_model10$adjR2
df_kastner_model$logLik_both <- df_modelBoth$logLik
df_kastner_model$AIC_both <- df_modelBoth$AIC
df_kastner_model$adjR2_both <- df_modelBoth$adjR2
df_kastner_model$both_ratio <- df_modelBoth$nonLin

# get roiname, side and subject properly
getRoiName <- function( idxLine, dataIn) { strsplit( as.character( dataIn$roiName[ idxLine ] ), '_' )[[1]][1] }
df_kastner_model$roi <- sapply( 1:dim(df_kastner_model)[1], getRoiName, dataIn=df_kastner_model )
getRoiSide <- function( idxLine, dataIn) { strsplit( as.character( dataIn$roiName[ idxLine ] ), '_' )[[1]][2] }
df_kastner_model$side <- sapply( 1:dim(df_kastner_model)[1], getRoiSide, dataIn=df_kastner_model )
getRoiSubject <- function( idxLine, dataIn) { strsplit( as.character( dataIn$subj[ idxLine ] ), '_' )[[1]][1] }
df_kastner_model$subject <- sapply( 1:dim(df_kastner_model)[1], getRoiSubject, dataIn=df_kastner_model )
df_kastner_model$roi <- ifelse( df_kastner_model$roi=='LO', 'LO', 
                                ifelse( df_kastner_model$roi=='LO1', 'LO', 
                                        ifelse( df_kastner_model$roi=='LO2', 'LO', df_kastner_model$roi ) ) )
df_kastner_model$roi <- ifelse( df_kastner_model$roi=='MT', 'MT', 
                                ifelse( df_kastner_model$roi=='MT1', 'MT', 
                                        ifelse( df_kastner_model$roi=='MT2', 'MT', df_kastner_model$roi ) ) )
df_kastner_model$roi <- ifelse( df_kastner_model$roi=='V2d', 'V2', 
                                ifelse( df_kastner_model$roi=='V2v', 'V2', df_kastner_model$roi ) )
df_kastner_model$roi <- ifelse( df_kastner_model$roi=='V3d', 'V3', 
                                ifelse( df_kastner_model$roi=='V3v', 'V3', df_kastner_model$roi ) )
df_kastner_model$roi_coarse <- ifelse( df_kastner_model$roi=='V1' | df_kastner_model$roi=='V2' | df_kastner_model$roi=='V3' | df_kastner_model$roi=='V4', 'V1-4', 
                                       ifelse( df_kastner_model$roi=='V3A' | df_kastner_model$roi=='V3B', 'V3AB', 
                                               ifelse( df_kastner_model$roi=='LO' | df_kastner_model$roi=='MT', 'LO/MT', 
                                                       ifelse( df_kastner_model$roi=='IPS0' | df_kastner_model$roi=='IPS1' | df_kastner_model$roi=='IPS2', 'IPS0-2', 
                                                               ifelse( df_kastner_model$roi=='IPS3' | df_kastner_model$roi=='IPS4' | df_kastner_model$roi=='IPS5', 'IPS3-5', 'FEF-PEF' ) ) ) ) )
table(df_kastner_model$roi_coarse, df_kastner_model$roi)

df_kastner_model$roi_coarse <- as.factor( df_kastner_model$roi_coarse )
df_kastner_model$roi_coarse <- factor( df_kastner_model$roi_coarse, levels=levels( df_kastner_model$roi_coarse )[c(5,6,4,2,3,1)] )

df_kastner_model_filt <- df_kastner_model[ df_prf$radius > 0.5 &
                                             df_prf$radius < 10 &
                                             df_prf$sigmaPos > 0 &
                                             df_prf$sigmaPos < 15 &
                                             df_prf$roi!='FEF' & df_prf$roi!='PEF', ]

cutData <- function( subjId, dataIn, roiInd, varThr, dataKastner ) {
  dataSel <- dataIn[ dataIn$subject==subjId & 
                       dataIn$roi==roiInd, ]
  dataSelKastner <- dataKastner[ dataKastner$subject==subjId & 
                                   dataKastner$roi==roiInd, ]
  varCut <- dataSel$var > quantile( dataSel$var, varThr ) # based on pRF data (dataSel is the prf dataset)
  keepIdx <- which( varCut )
  dataSel <- dataSelKastner[ varCut, ]
  dataSel$keepIdx <- keepIdx
  return( dataSel )
}

# for each participant and roi, cut the data
flagDf <- 1
for ( participantId in unique( df_kastner_model_filt$subject ) ) {
  df_filt_subj <- df_kastner_model_filt[ df_kastner_model_filt$subject==participantId, ]
  df_filt_subj_prf <- df_prf_filt[ df_prf_filt$subject==participantId, ]
  for ( roiId in unique( df_filt_subj$roi ) ) {
    dataDfTemp <- cutData( participantId, df_prf_filt, roiId, 0.50, df_kastner_model_filt  ) #here is the line to filter the dataset based on pRF data
    if (flagDf==1) { df_kastner_model_filt_clean <- dataDfTemp; flagDf <- 0  }
    if (flagDf==0) { df_kastner_model_filt_clean <- rbind( df_kastner_model_filt_clean, dataDfTemp ) }
  }
}
str( df_kastner_model_filt_clean ) 
table( df_kastner_model_filt_clean$roi, df_kastner_model_filt_clean$subject )

# var explained
str( df_kastner_model )
df_kastner_model_aggregate_var <- aggregate( varExp ~ subject * roi, mean, data= df_kastner_model_filt_clean ) 
roiLevels <- unique(  df_kastner_model_aggregate_var$roi )
roiLevels <- roiLevels[ c(9,10,11,12,13,14,7,8,1,2,3,4,5,6) ] 
df_kastner_model_aggregate_var$roi_sorted <- factor( df_kastner_model_aggregate_var$roi, levels=roiLevels ) 
boxplot( df_kastner_model_aggregate_var$varExp ~ df_kastner_model_aggregate_var$roi_sorted, ylim=c(0,0.4), frame=FALSE, las=2, col=subjCols, outline=FALSE )
summary( lm( varExp ~ roi_sorted, data=df_kastner_model_aggregate_var) )

# kappa par
df_kastner_model_aggregate_kappa <- aggregate( kappa ~ subject * roi, mean, data= df_kastner_model_filt_clean ) 
roiLevels <- unique(  df_kastner_model_aggregate_kappa$roi )
roiLevels <- roiLevels[ c(9,10,11,12,13,14,7,8,1,2,3,4,5,6) ] 
df_kastner_model_aggregate_kappa$roi_sorted <- factor( df_kastner_model_aggregate_kappa$roi, levels=roiLevels ) 
boxplot( df_kastner_model_aggregate_kappa$kappa ~ df_kastner_model_aggregate_kappa$roi_sorted, ylim=c(0.2,0.8), frame=FALSE, las=2, col=subjCols, outline=FALSE  )
summary( lm( kappa ~ roi_sorted, data=df_kastner_model_aggregate_kappa) )

# all rois, adjR2 diff wrt model 10
df_kastner_model_aggregate_adjR2DF <- aggregate( adjR2 ~ subject * roi, mean, data= df_kastner_model_filt_clean ) 
df_kastner_model_aggregate_adjR2DF10 <- aggregate( adjR2_10 ~ subject * roi, mean, data= df_kastner_model_filt_clean ) 
df_kastner_model_aggregate_adjR2DF$adjR2_10 <- df_kastner_model_aggregate_adjR2DF10$adjR2_10
df_kastner_model_aggregate_adjR2DF$adjR2_diff <- df_kastner_model_aggregate_adjR2DF$adjR2 - df_kastner_model_aggregate_adjR2DF10$adjR2_10
roiLevels <- unique(  df_kastner_model_aggregate_adjR2DF$roi )
roiLevels <- roiLevels[ c(9,10,11,12,13,14,7,8,1,2,3,4,5,6) ] 
df_kastner_model_aggregate_adjR2DF$roi_sorted <- factor( df_kastner_model_aggregate_adjR2DF$roi, levels=roiLevels ) 
par(mfrow=c(1,1))
boxplot( df_kastner_model_aggregate_adjR2DF$adjR2_diff ~ df_kastner_model_aggregate_adjR2DF$roi_sorted, frame=FALSE, las=1, col=subjCols, ylim=c(-0.1,0.4)  )
summary( lm( adjR2_diff ~ roi_sorted, data=df_kastner_model_aggregate_adjR2DF) )

# all rois, adjR2 diff wrt model both
df_kastner_model_aggregate_adjR2DF <- aggregate( adjR2 ~ subject * roi, mean, data= df_kastner_model_filt_clean ) 
df_kastner_model_aggregate_adjR2DF10 <- aggregate( adjR2_both ~ subject * roi, mean, data= df_kastner_model_filt_clean ) 
df_kastner_model_aggregate_adjR2DF$adjR2_10 <- df_kastner_model_aggregate_adjR2DF10$adjR2_both
df_kastner_model_aggregate_adjR2DF$adjR2_diff <- df_kastner_model_aggregate_adjR2DF$adjR2 - df_kastner_model_aggregate_adjR2DF10$adjR2_both
roiLevels <- unique(  df_kastner_model_aggregate_adjR2DF$roi )
roiLevels <- roiLevels[ c(9,10,11,12,13,14,7,8,1,2,3,4,5,6) ] 
df_kastner_model_aggregate_adjR2DF$roi_sorted <- factor( df_kastner_model_aggregate_adjR2DF$roi, levels=roiLevels ) 
par(mfrow=c(1,1))
boxplot( df_kastner_model_aggregate_adjR2DF$adjR2_diff ~ df_kastner_model_aggregate_adjR2DF$roi_sorted, frame=FALSE, las=1, col=subjCols, ylim=c(-0.1,0.4)  )
summary( lm( adjR2_diff ~ roi_sorted, data=df_kastner_model_aggregate_adjR2DF) )

# all rois, both multiplicative par
df_kastner_model_aggregate_both_ratio <- aggregate( both_ratio ~ subject * roi, mean, data= df_kastner_model_filt_clean ) 
roiLevels <- unique(  df_kastner_model_aggregate_both_ratio$roi )
roiLevels <- roiLevels[ c(9,10,11,12,13,14,7,8,1,2,3,4,5,6) ] 
df_kastner_model_aggregate_both_ratio$roi_sorted <- factor( df_kastner_model_aggregate_both_ratio$roi, levels=roiLevels ) 
par(mfrow=c(1,1))
boxplot( df_kastner_model_aggregate_both_ratio$both_ratio ~ df_kastner_model_aggregate_both_ratio$roi_sorted, ylim=c(0,1), frame=FALSE, las=1, col=subjCols  )
summary( lm( both_ratio ~ roi_sorted, data=df_kastner_model_aggregate_both_ratio) )

# with histograms, density
par( mfrow=c(1,5) )
str( df_kastner_model_filt_clean )
nBreaks <- 25
dfFlag <- 1
for ( subjSel in unique( df_kastner_model_filt_clean$subject )  ) {
  for ( roiSel in levels( df_kastner_model_filt_clean$roi_coarse )[1:5]  ) {
    
    df_kastner_model_filt_subj_left <- subset( df_kastner_model_filt_clean, subject==subjSel & roi_coarse==roiSel & side=='left' )
    ds_left <- hist( df_kastner_model_filt_subj_left$mu, plot=FALSE, breaks = seq(0,2*pi, length.out=nBreaks) )
    df_kastner_model_filt_subj_right <- subset( df_kastner_model_filt_clean, subject==subjSel & roi_coarse==roiSel & side=='right' )
    ds_right <- hist( df_kastner_model_filt_subj_right$mu, plot=FALSE, breaks = seq(0,2*pi, length.out=nBreaks) )
    
    # density
    maxD <- max( c( max(ds_left$density), max(ds_right$density) )  )
    plot( 0, type='n', ylab='density', xlab='mu (rad)', main=roiSel, xlim=c(0,2*pi), ylim=c(0,maxD), axes=FALSE, bty='n' )
    lines( ds_left$mids, ds_left$density , col='darkorange', lwd=3, lty=1 )
    lines( ds_right$mids, ds_right$density, col='darkblue', lwd=3, lty=1 )
    axis( 1, c(0,pi,2*pi), round( c(0,pi,2*pi),2 ) ); axis(2, c(0,round(maxD,2)), las=1 )
    
    subjArrayLeft <-  rep( subjSel, length( ds_left$density ) )
    sideArrayLeft <- rep( 'left', length( ds_left$density ) )
    roiArrayLeft <-  rep( roiSel, length( ds_left$density ) )
    densityArrayLeft <-  ds_left$density
    midsArrayLeft <- ds_left$mids
    
    subjArrayRight <-  rep( subjSel, length( ds_right$density ) )
    sideArrayRight <- rep( 'right', length( ds_right$density ) )
    roiArrayRight <-  rep( roiSel, length( ds_right$density ) )
    densityArrayRight <-  ds_right$density
    midsArrayRight <- ds_right$mids
    
    subjArray <- c( subjArrayLeft, subjArrayRight )
    sideArray <- c( sideArrayLeft, sideArrayRight )
    roiArray <- c( roiArrayLeft, roiArrayRight )
    densityArray <- c( densityArrayLeft, densityArrayRight )
    midsArray <- c( midsArrayLeft, midsArrayRight )
    
    subjDF <- data.frame( subjArray, sideArray, roiArray, densityArray, midsArray  )
    
    if (dfFlag == 1) { subjDFOutDensity <- subjDF; dfFlag <- 2 }
    if (dfFlag>1) { subjDFOutDensity <- rbind( subjDFOutDensity, subjDF ) } 
    
  }
}

aggDFDensity <- aggregate( densityArray ~ midsArray*sideArray*roiArray, mean, data=subjDFOutDensity )
aggDFDensitySd <- aggregate( densityArray ~ midsArray*sideArray*roiArray, sd, data=subjDFOutDensity )
aggDFDensitySd$Se <- aggDFDensitySd$densityArray  / sqrt( length( unique( df_kastner_model_filt_clean$subject ) ) )

par( mfrow=c(1,5) )
for ( roiIdx in levels( aggDFDensity$roiArray ) ) {
  selectedDf <- subset( aggDFDensity, roiArray==roiIdx )
  selectedDfSd <- subset( aggDFDensitySd, roiArray==roiIdx )
  maxD <- max( c( max(selectedDf$densityArray), max(selectedDf$densityArray) )  )
  maxD <- maxD + maxD*0.15
  plot( 0, type='n', ylab='density', xlab='mu (rad)', main=roiIdx, xlim=c(0,2*pi), ylim=c(0,maxD), axes=FALSE, bty='n' )
  lines( selectedDf$midsArray[ selectedDf$sideArray=='left' ], selectedDf$densityArray[ selectedDf$sideArray=='left' ], col='darkorange', lwd=3, lty=1 )
  lines( selectedDf$midsArray[ selectedDf$sideArray=='left' ], selectedDf$densityArray[ selectedDf$sideArray=='left' ]+selectedDfSd$Se[ selectedDf$sideArray=='left' ], col='orange', lwd=1, lty=1 )
  lines( selectedDf$midsArray[ selectedDf$sideArray=='left' ], selectedDf$densityArray[ selectedDf$sideArray=='left' ]-selectedDfSd$Se[ selectedDf$sideArray=='left' ], col='orange', lwd=1, lty=1 )
  
  lines( selectedDf$midsArray[ selectedDf$sideArray=='right' ], selectedDf$densityArray[ selectedDf$sideArray=='right' ], col='darkblue', lwd=3, lty=1 )
  lines( selectedDf$midsArray[ selectedDf$sideArray=='right' ], selectedDf$densityArray[ selectedDf$sideArray=='right' ]+selectedDfSd$Se[ selectedDf$sideArray=='right' ], col='darkblue', lwd=1, lty=1 )
  lines( selectedDf$midsArray[ selectedDf$sideArray=='right' ], selectedDf$densityArray[ selectedDf$sideArray=='right' ]-selectedDfSd$Se[ selectedDf$sideArray=='right' ], col='darkblue', lwd=1, lty=1 )
  axis( 1, c(0,pi,2*pi), round( c(0,pi,2*pi),2 ) ); axis(2, c(0,round(maxD,2)), las=1 )
}


#### add the analysis for phase encoded design, with both thresholds
