rm( list=ls() );
mainDir <- '/analyse/Project0226/KastnerModel/staging_area_Kastner'
pRFDir <- 'pRF'
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

#### load data, model pRF ####
load( sprintf('%s/_pRF_params.nii.gz.RData', pRFDir ) )
df_model_pRF <- outModelDf
rm( list = c( 'outModelDf' ) )
str( df_model_pRF )

df_model_pRF$radius
df_model_pRF$participant
df_model_pRF$ROI_name_basic
df_model_pRF$varExp

# cut data function
dataIn <- df_model_pRF
participantId <- 'AFH28'
propCut <- c(0,1,0.1)
roiInd <- 'V1'
varThr <- 0.6

cutData <- function( dataIn, participantId, propCut, roiInd, varThr ) {
  dataSel <- dataIn[ dataIn$participant==participantId & 
                       dataIn$ROI_name_basic==roiInd, ]
  varCut <- dataSel$varExp > quantile( dataSel$varExp, varThr )
  #varCut <- dataSel$varExp > varThr
  
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
    dataDfTemp <- cutData(df_model_pRF, participantId, c(0,1,0.05), roiId, 0.65  ) 
    if (flagDf==1) { df_prf_filt_clean_full <- dataDfTemp; flagDf <- 0  }
    if (flagDf==0) { df_prf_filt_clean_full <- rbind( df_prf_filt_clean_full, dataDfTemp ) }
  }
}
str( df_prf_filt_clean_full ) 
table( df_prf_filt_clean_full$ROI_name_basic, df_prf_filt_clean_full$participant )

##### select only ROIs of interest #####

df_prf_filt_clean_full$ROI_name <- as.character( df_prf_filt_clean_full$ROI_name )
#df_prf_filt_clean <- subset( df_prf_filt_clean_full, ROI_name=='V1v' | ROI_name=='V1d' | ROI_name=='V2v'  | ROI_name=='V2d' 
#                             | ROI_name=='V3v' | ROI_name=='V3d' )
df_prf_filt_clean <- subset( df_prf_filt_clean_full, ROI_name=='V1v' | ROI_name=='V1d' )
roiLevels <- unique(  df_prf_filt_clean$ROI_name )
df_prf_filt_clean$roi_sorted <- factor( df_prf_filt_clean$ROI_name, levels=roiLevels ) 
df_prf_filt_clean$ROI_name <- df_prf_filt_clean$roi_sorted

df_prf_filt_clean$dv <- df_prf_filt_clean$sigmaPos
aggData <- aggregate( dv ~ roi_sorted * participant, median, data=df_prf_filt_clean )
boxplot( aggData$dv ~ aggData$roi_sorted, notch=TRUE )
barplot( tapply( aggData$dv, list( aggData$roi_sorted ), mean ) )
tapply( aggData$dv, list( aggData$roi_sorted ), mean )
t.test( aggData$dv[aggData$roi_sorted=='V1v'], aggData$dv[aggData$roi_sorted=='V1d'], paired=TRUE)

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

  plot( sigmaPos ~ radius, data=size_ecc_df_base[ size_ecc_df_base$roi_sorted=='V1v' & size_ecc_df_base$participant==participantId, ], 
        col='red', xlab='eccentricity (dva)', ylab='size (dva)', cex.lab=1.55, cex.axis=1.5,
        ylim=c(0,4), xlim=c(0,10), pch=8, las=1, bty='n' ); 
  abline( lm( sigmaPos ~ radius, data=size_ecc_df_base[ size_ecc_df_base$roi_sorted=='V1v' & size_ecc_df_base$participant==participantId, ] ), col='red' )
  text( 2.5, 9, participantId, cex=1.25 )
  
  dataTemp <- size_ecc_df_base[ size_ecc_df_base$roi_sorted=='V1d' & size_ecc_df_base$participant==participantId, ]
  points( dataTemp$sigmaPos~dataTemp$radius, col='black');
  abline( lm( sigmaPos ~ radius, data=size_ecc_df_base[ size_ecc_df_base$roi_sorted=='V1d' & size_ecc_df_base$participant==participantId, ] ), col='black' )
  
}

# size x eccentricity, global average, plot
plot( sigmaPos ~ radius, data=size_ecc_df_base_globalAverage[ size_ecc_df_base_globalAverage$roi_sorted=='V1v', ], col='red', 
      ylim=c(0,3), xlim=c(0,10), xlab='eccentricity (dva)', ylab='size (dva)', cex.lab=1.55, cex.axis=1.5,
      pch=8, las=1, bty='n' ); abline( lm( sigmaPos ~ radius, data=size_ecc_df_base_globalAverage[ size_ecc_df_base_globalAverage$roi_sorted=='V1v', ] ), col='red' )
text( 3.2, 9, 'Average', cex=1.5 )
dataTemp <- size_ecc_df_base_globalAverage[ size_ecc_df_base_globalAverage$roi_sorted=='V1d', ]
points( dataTemp$sigmaPos~dataTemp$radius, col='black');
abline( lm( sigmaPos ~ radius, data=size_ecc_df_base_globalAverage[ size_ecc_df_base$roi_sorted=='V1d', ] ), col='black' )

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
size_ecc_df_base_globalAverage_var <- aggregate( varExp ~ participant * roi_sorted, mean, data=size_ecc_df_base)  #here varExp refers to adjR2, see line 88 and line 90
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
