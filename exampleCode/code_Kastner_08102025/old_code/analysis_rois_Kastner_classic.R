rm( list=ls() );
resultsDir <- '/analyse/Project0226/KastnerModel/results_KastnerClassic'
filenameOut <- 'testResults_all_in_Kastner_classic_withoutHrf_30092021.RData'
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

#### load data ####
load( filenameOut )
str( df_prf )
str( df_phase10_CCW )

#####  get roiname, side and subject properly, prepare data and declare cut data function ##### 
getRoiName <- function( idxLine, dataIn) { strsplit( as.character( dataIn$roiName[ idxLine ] ), '_' )[[1]][1] }
df_prf$roi <- sapply( 1:dim(df_prf)[1], getRoiName, dataIn=df_prf )
getRoiSide <- function( idxLine, dataIn) { strsplit( as.character( dataIn$roiName[ idxLine ] ), '_' )[[1]][2] }
df_prf$side <- sapply( 1:dim(df_prf)[1], getRoiSide, dataIn=df_prf )
getRoiSubject <- function( idxLine, dataIn) { strsplit( as.character( dataIn$subj[ idxLine ] ), '_' )[[1]][1] }
df_prf$subject <- sapply( 1:dim(df_prf)[1], getRoiSubject, dataIn=df_prf )
df_prf$roi <- ifelse( df_prf$roi=='LO', 'LO', 
                           ifelse( df_prf$roi=='LO1', 'LO', 
                                   ifelse( df_prf$roi=='LO2', 'LO', df_prf$roi ) ) )
df_prf$roi <- ifelse( df_prf$roi=='MT', 'MT', 
                      ifelse( df_prf$roi=='MT1', 'MT', 
                              ifelse( df_prf$roi=='MT2', 'MT', df_prf$roi ) ) )
df_prf$roi <- ifelse( df_prf$roi=='V2d', 'V2', 
                      ifelse( df_prf$roi=='V2v', 'V2', df_prf$roi ) )
df_prf$roi <- ifelse( df_prf$roi=='V3d', 'V3', 
                      ifelse( df_prf$roi=='V3v', 'V3', df_prf$roi ) )

df_prf_filt <- df_prf[ df_prf$radius > 0.5 &
                         df_prf$radius < 10 &
                         df_prf$sigmaPos > 0 &
                         df_prf$sigmaPos < 15 &
                         df_prf$roi!='FEF' & df_prf$roi!='PEF', ]

table( df_prf_filt$roi, df_prf_filt$subject )

# adjusted r squared for prf data
nParams_prf <- 3
df_prf_filt$adjR2 <- 1 - (1-df_prf_filt$varExp) * ( (155-1) / (155-nParams_prf-1) )

# cut data function
cutData <- function( subjId, dataIn, propCut, roiInd, varThr ) {
  dataSel <- dataIn[ dataIn$subject==subjId & 
                       dataIn$roi==roiInd, ]
  varCut <- dataSel$var > quantile( dataSel$adjR2, varThr )
  #varCut <- dataSel$var > varThr
  keptIdx <- which( varCut )
  dataSel <- dataSel[ varCut, ]
  eccCut <- .bincode( dataSel$radius, 
                      breaks = quantile( dataSel$radius, probs =  seq( propCut[1], propCut[2], propCut[3] ) ), 
                      include.lowest = TRUE )
  xCut <- .bincode( dataSel$xPos, 
                      breaks = quantile( dataSel$xPos, probs =  seq( propCut[1], propCut[2], propCut[3] ) ), 
                      include.lowest = TRUE )
  yCut <- .bincode( dataSel$yPos, 
                    breaks = quantile( dataSel$yPos, probs =  seq( propCut[1], propCut[2], propCut[3] ) ), 
                    include.lowest = TRUE )
  dataSel$eccCut <- as.numeric( eccCut )
  dataSel$xCut <- as.numeric( xCut )
  dataSel$yCut <- as.numeric( yCut )
  dataSel$keptIdx <- keptIdx
  return( dataSel )
}

##### for each participant and roi, cut the data #####
flagDf <- 1
for ( participantId in unique( df_prf_filt$subject ) ) {
  df_prf_filt_subj <- df_prf_filt[ df_prf_filt$subject==participantId, ]
  for ( roiId in unique( df_prf_filt_subj$roi ) ) {
    dataDfTemp <- cutData( participantId, df_prf_filt, c(0,1,0.05), roiId, 0.65  ) 
    if (flagDf==1) { df_prf_filt_clean <- dataDfTemp; flagDf <- 0  }
    if (flagDf==0) { df_prf_filt_clean <- rbind( df_prf_filt_clean, dataDfTemp ) }
  }
}
str( df_prf_filt_clean ) 
table( df_prf_filt_clean$roi, df_prf_filt_clean$subject )

##### size x eccentricy, single subjects, plot and analysis ##### 

# size x eccentricy, single subjects
size_ecc_df_base <- aggregate( sigmaPos ~ eccCut * roi * subject, mean, data=df_prf_filt_clean )
size_ecc_df_radius <- aggregate( radius ~ eccCut * roi * subject, mean, data=df_prf_filt_clean )
size_ecc_df_varExp <- aggregate( adjR2 ~ eccCut * roi * subject, mean, data=df_prf_filt_clean )
size_ecc_df_base$radius <- size_ecc_df_radius$radius
size_ecc_df_base$varExp <- size_ecc_df_varExp$adjR2

# size x eccentricity, global average
size_ecc_df_base_globalAverage <- aggregate( sigmaPos ~ eccCut * roi, mean, data=size_ecc_df_base) 
size_ecc_df_radius_globalAverage <- aggregate( radius ~ eccCut * roi, mean, data=size_ecc_df_base) 
size_ecc_df_base_globalAverage$radius <- size_ecc_df_radius_globalAverage$radius

x11( width=7, height=8.1 )

par(mfrow=c(3,3))
subjCols <- rainbow( length( unique( df_prf_filt_clean$roi ) ), star=0, end=0.80 )
for ( participantId in unique( df_prf_filt_clean$subject ) ) {
  plot( sigmaPos ~ radius, data=size_ecc_df_base[ size_ecc_df_base$roi=='V1' & size_ecc_df_base$subject==participantId, ], 
        col=subjCols[ 1 ], xlab='eccentricity (dva)', ylab='size (dva)', cex.lab=1.55, cex.axis=1.5,
        ylim=c(0,10), xlim=c(0,10), pch=8, las=1, bty='n' ); 
  abline( lm( sigmaPos ~ radius, data=size_ecc_df_base[ size_ecc_df_base$roi=='V1' & size_ecc_df_base$subject==participantId, ] ), col=subjCols[ 1 ] )
  text( 2, 9, participantId )
  roiLevels <- unique( size_ecc_df_base$roi )
  roiLevels <- roiLevels[ c(9,10,11,12,13,14,7,8,1,2,3,4,5,6) ] 
  roiCounter <- 2
  for (roiIdx in roiLevels ) {
    dataPlot <- size_ecc_df_base[ size_ecc_df_base$roi==roiIdx & size_ecc_df_base$subject==participantId, ]; points( dataPlot$radius, dataPlot$sigmaPos, col=subjCols[ roiCounter ], pch=8 ); abline( lm( dataPlot$sigmaPos ~ dataPlot$radius ), col=subjCols[ roiCounter ] )
    roiCounter <- roiCounter + 1
  }
}

# size x eccentricity, global average, plot
plot( sigmaPos ~ radius, data=size_ecc_df_base_globalAverage[ size_ecc_df_base_globalAverage$roi=='V1', ], col=subjCols[ 1 ], 
      ylim=c(0,10), xlim=c(0,10), xlab='eccentricity (dva)', ylab='size (dva)', cex.lab=1.55, cex.axis=1.5,
      pch=8, las=1, bty='n' ); abline( lm( sigmaPos ~ radius, data=size_ecc_df_base_globalAverage[ size_ecc_df_base_globalAverage$roi=='V1', ] ), col=subjCols[ 1 ] )
text( 2.5, 9, 'Average' )
roiLevels <- unique( size_ecc_df_base_globalAverage$roi )
roiLevels <- roiLevels[ c(9,10,11,12,13,14,7,8,1,2,3,4,5,6) ] 
roiCounter <- 2
for (roiIdx in roiLevels ) {
  dataPlot <- size_ecc_df_base_globalAverage[ size_ecc_df_base_globalAverage$roi==roiIdx, ]; points( dataPlot$radius, dataPlot$sigmaPos, col=subjCols[ roiCounter ], pch=8 ); abline( lm( dataPlot$sigmaPos ~ dataPlot$radius ), col=subjCols[ roiCounter ] )
  roiCounter <- roiCounter + 1
}

#legend plot
plot(0,axes=FALSE,bty='n',xlab='',ylab='',type='n',xlim=c(0,1),ylim=c(0,1))
legend(0, 1.0,legend=roiLevels[1:5],fill=subjCols[1:5],bty='n',cex=0.8)
legend(0.4, 1.0,legend=roiLevels[6:10],fill=subjCols[6:10],bty='n',cex=0.8)
legend(0.7, 1.0,legend=roiLevels[11:14],fill=subjCols[11:14],bty='n',cex=0.8)

dev.copy2pdf( file='prf_sizeEccentricity.pdf', width=7, height=7 )
graphics.off()

size_ecc_df_base$roi_factor <- as.factor( size_ecc_df_base$roi )
size_ecc_df_base$eccCut_scalar <- as.numeric( size_ecc_df_base$eccCut )
size_ecc_df_base$subject_factor <- as.factor( size_ecc_df_base$subject )
summary( lmer( sigmaPos ~ eccCut_scalar + roi_factor + ( 1 | subject_factor ), data=size_ecc_df_base ) )

##### prf variance explained plot and analysis ##### 
x11(width=4.5, height=4.5)

par(mfrow=c(1,1))
size_ecc_df_base_globalAverage_var <- aggregate( varExp ~ subject * roi, mean, data=size_ecc_df_base)  #here varExp refers to adjR2, see line 88 and line 90
roiLevels <- unique( size_ecc_df_base_globalAverage_var$roi )
roiLevels <- roiLevels[ c(9,10,11,12,13,14,7,8,1,2,3,4,5,6) ] 
size_ecc_df_base_globalAverage_var$roi_sorted <- factor( size_ecc_df_base_globalAverage_var$roi, levels=roiLevels ) 
boxplot( size_ecc_df_base_globalAverage_var$varExp ~ size_ecc_df_base_globalAverage_var$roi_sorted, ylim=c(0,1), 
         frame=FALSE, las=2, col=subjCols, outline=FALSE, xlab='ROI', ylab='var. exp. (proportion)', cex.lab=1.35, cex.axis=1.25 )

dev.copy2pdf( file='prf_varExp.pdf', width=4.5, height=4.5 )
graphics.off()

summary( lmer( varExp ~ roi_sorted + ( 1 | subject ), data=size_ecc_df_base_globalAverage_var ) )


# # up VS down analysis x eccentricity
# str(df_prf_filt_clean)
# df_prf_filt_clean$up_down <- ifelse( df_prf_filt_clean$yPos<0,'LVF','UVF')
# df_xyposition <- aggregate( sigmaPos ~ subj * roi * eccCut * up_down, FUN=mean, data=df_prf_filt_clean )
# df_xypositionXpos <- aggregate( xPos ~ subj * roi * eccCut * up_down, FUN=mean, data=df_prf_filt_clean )
# df_xypositionYpos <- aggregate( yPos ~ subj * roi * eccCut * up_down, FUN=mean, data=df_prf_filt_clean )
# df_xypositionEcc <- aggregate( radius ~ subj * roi * eccCut * up_down, FUN=mean, data=df_prf_filt_clean )
# df_xyposition$xPos <- df_xypositionXpos$xPos
# df_xyposition$yPos <- df_xypositionYpos$yPos
# df_xyposition$ecc <- df_xypositionEcc$radius
# df_xyposition_aggregate <- aggregate( sigmaPos ~ roi * eccCut * up_down, FUN=mean, data=df_xyposition )
# df_xyposition_aggregateEcc <- aggregate( ecc ~ roi * eccCut * up_down, FUN=mean, data=df_xyposition )
# df_xyposition_aggregate$ecc <- df_xyposition_aggregateEcc$ecc 
# # size x eccentricity, global average, plot
# par(mfrow=c(1,1))
# plot( sigmaPos ~ ecc, data=df_xyposition_aggregate[ df_xyposition_aggregate$roi=='V1' & df_xyposition_aggregate$up_down=='LVF', ], col='darkblue', ylim=c(0,3), xlim=c(0,10), pch=8, las=1, bty='n' ); 
# abline( lm( sigmaPos ~ radius, data=size_ecc_df_base_globalAverage[ size_ecc_df_base_globalAverage$roi=='V1' & df_xyposition_aggregate$up_down=='LVF', ] ), col='darkblue' )
# points( df_xyposition_aggregate$ecc[ df_xyposition_aggregate$roi=='V1' & df_xyposition_aggregate$up_down=='UVF' ], df_xyposition_aggregate$sigmaPos[ df_xyposition_aggregate$roi=='V1' & df_xyposition_aggregate$up_down=='UVF' ], col='darkred', pch=8 ); 
# abline( lm( df_xyposition_aggregate$sigmaPos[ df_xyposition_aggregate$roi=='V1' & df_xyposition_aggregate$up_down=='UVF' ] ~ df_xyposition_aggregate$ecc[ df_xyposition_aggregate$roi=='V1' & df_xyposition_aggregate$up_down=='UVF' ] ), col='darkred' )
# 
# # up VS down analysis, simple, V1
# str(df_prf_filt_clean)
# df_prf_filt_clean$up_down <- ifelse( df_prf_filt_clean$yPos<0,'LVF','UVF')
# df_xyposition <- aggregate( sigmaPos ~ subj * roi * up_down, FUN=mean, data=df_prf_filt_clean )
# par(mfrow=c(1,2))
# boxplot( df_xyposition$sigmaPos[ df_xyposition$roi=='V1' & df_xyposition$up_down=='LVF' ],
#          df_xyposition$sigmaPos[ df_xyposition$roi=='V1' & df_xyposition$up_down=='UVF' ], names=c('LVF','UVF'), frame=FALSE, las=1, outline = FALSE, ylim=c(1.2,2))
# yLVF <- df_xyposition$sigmaPos[ df_xyposition$roi=='V1' & df_xyposition$up_down=='LVF' ]
# yUVF <- df_xyposition$sigmaPos[ df_xyposition$roi=='V1' & df_xyposition$up_down=='UVF' ]
# plot( 0, xlim=c(0,2.5), ylim=c(0.2,2.5), bty='n', las=1, ylab='pRF size (dva)', xlab='', axes=FALSE )
# points( rep(0.5,length(yLVF)), yLVF, cex=2, pch=seq(8,8+length(yLVF)) ); points( rep(2,length(yLVF)), yUVF, cex=2, pch=seq(8,8+length(yLVF)) )
# segments( rep(0.5,length(yLVF)), yLVF, rep(2,length(yUVF)), yUVF, lwd=2, lty=2 )
# axis( 2, seq(0.2,2.5,length.out=3 ), las=1 ); axis( 1, c(0.5,2), c('LVF','UVF') )
# t.test(  df_xyposition$sigmaPos[ df_xyposition$roi=='V1' & df_xyposition$up_down=='LVF' ],
#          df_xyposition$sigmaPos[ df_xyposition$roi=='V1' & df_xyposition$up_down=='UVF' ], paired=TRUE )
# 
# # up VS down analysis, xCut, yCut
# par(mfrow=c(1,1))
# df_xyposition <- aggregate( sigmaPos ~ subj * roi * xCut * yCut, FUN=mean, data=df_prf_filt_clean )
# df_xypositionXpos <- aggregate( xPos ~ subj * roi * xCut * yCut, FUN=mean, data=df_prf_filt_clean )
# df_xypositionYpos <- aggregate( yPos ~ subj * roi * xCut * yCut, FUN=mean, data=df_prf_filt_clean )
# df_xyposition$xPos <- df_xypositionXpos$xPos
# df_xyposition$yPos <- df_xypositionYpos$yPos
# plot( df_xyposition$yPos ~ df_xyposition$xPos, bty='n', las=1 )
# df_xyposition$up_down <- ifelse( df_xyposition$yPos<0,'LVF','UVF')
# str( df_xyposition )
# df_updown <- aggregate( sigmaPos ~ subj * roi * up_down, FUN=mean, data=df_xyposition )
# df_updown$sigmaPos[ df_updown$roi=='V1' & df_updown$up_down=='LVF' ]
# df_updown$sigmaPos[ df_updown$roi=='V1' & df_updown$up_down=='UVF' ]
# t.test( df_updown$sigmaPos[ df_updown$roi=='V1' & df_updown$up_down=='LVF' ],
#         df_updown$sigmaPos[ df_updown$roi=='V1' & df_updown$up_down=='UVF' ], paired=TRUE )



##################################################################
##################################################################
## model kastner (filtered on model_kastner classic data itself)##
##################################################################
##################################################################
str(df_model05)
df_kastner_model <- df_model05; 

df_kastner_model$logLik_05 <- df_kastner_model$logLik; nParams_05 <- 4 # location / width / intercept / slope
df_kastner_model$AIC_05 <- 2*nParams_05 - 2*df_kastner_model$logLik_05
df_kastner_model$adjR2_05 <- 1 - (1-df_kastner_model$varExp) * ( (279-1) / (279-nParams_05-1) )
df_kastner_model$mu05_CW <- df_model05_CW$mu
df_kastner_model$mu05_CCW <- df_model05_CCW$mu
df_kastner_model$kappa05_CW <- df_model05_CW$kappa
df_kastner_model$kappa05_CCW <- df_model05_CCW$kappa

df_kastner_model$logLik_10 <- df_model10$logLik; nParams_10 <- 4 # location / width / intercept / slope
df_kastner_model$AIC_10 <- 2*nParams_10 - 2*df_kastner_model$logLik_10 
df_kastner_model$adjR2_10 <- 1 - (1-df_model10$varExp) * ( (279-1) / (279-nParams_10-1) )
df_kastner_model$mu10_CW <- df_model10_CW$mu
df_kastner_model$mu10_CCW <- df_model10_CCW$mu
df_kastner_model$kappa10_CW <- df_model10_CW$kappa
df_kastner_model$kappa10_CCW <- df_model10_CCW$kappa

df_kastner_model$logLik_both <- df_modelBoth$logLik; nParams_both <- 5 # location / width / intercept / slope / ratio
df_kastner_model$AIC_both <- 2*nParams_both - 2*df_kastner_model$logLik_both
df_kastner_model$adjR2_both <- 1 - (1-df_modelBoth$varExp) * ( (279-1) / (279-nParams_both-1) )
df_kastner_model$both_ratio <- df_modelBoth$nonLin
df_kastner_model$muBoth_CW <- df_modelBoth_CW$mu
df_kastner_model$muBoth_CCW <- df_modelBoth_CCW$mu
df_kastner_model$kappaBoth_CW <- df_modelBoth_CCW$kappa
df_kastner_model$kappaBoth_CCW <- df_modelBoth_CCW$kappa

nParams_phase05 <- 2 # phase / amplitude 
nParams_phase10 <- 2 # phase / amplitude 
df_kastner_model$adjR2_phase05 <- 1 - ( 1-df_phase05$co^2 ) * ( (279-1) / (279-nParams_phase05-1) )
df_kastner_model$adjR2_phase10 <- 1 - ( 1-df_phase10$co^2 ) * ( (279-1) / (279-nParams_phase10-1) )
df_kastner_model$phase05 <- df_phase05$phase
df_kastner_model$phase05_CW <- ( df_phase05_CW$phase - 0.62 ) %% ( 2*pi )
df_kastner_model$phase05_CCW <- ( df_phase05_CCW$phase - 0.62 ) %% ( 2*pi )
df_kastner_model$phase10 <- df_phase10$phase
df_kastner_model$phase10_CW <- df_phase10_CW$phase
df_kastner_model$phase10_CCW <- df_phase10_CCW$phase
##
# from here, fix the phase problems with phase-encoded design, take weighted circular average from phase vectors
# perform roi analysis based on prf variance expained to avoid circularity
##

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

df_kastner_model_filt <- df_kastner_model[ df_kastner_model$roi!='FEF' & 
                                             df_kastner_model$roi!='PEF', ]

cutData <- function( subjId, dataIn, roiInd, varThr ) {
  dataSel <- dataIn[ dataIn$subject==subjId & 
                       dataIn$roi==roiInd, ]
  varCut <- dataSel$var > quantile( dataSel$var, varThr )
  dataSel <- dataSel[ varCut, ]
  return( dataSel )
}

# for each participant and roi, cut the data
flagDf <- 1
for ( participantId in unique( df_kastner_model_filt$subject ) ) {
  df_filt_subj <- df_kastner_model_filt[ df_kastner_model_filt$subject==participantId, ]
  for ( roiId in unique( df_filt_subj$roi ) ) {
    dataDfTemp <- cutData( participantId, df_kastner_model_filt, roiId, 0.65  ) 
    if (flagDf==1) { df_kastner_model_filt_clean <- dataDfTemp; flagDf <- 0  }
    if (flagDf==0) { df_kastner_model_filt_clean <- rbind( df_kastner_model_filt_clean, dataDfTemp ) }
  }
}
str( df_kastner_model_filt_clean ) 
table( df_kastner_model_filt_clean$roi, df_kastner_model_filt_clean$subject )

# all rois, var explained
str( df_kastner_model )
df_kastner_model_aggregate_var <- aggregate( adjR2_05 ~ subject * roi, mean, data= df_kastner_model_filt_clean ) 
roiLevels <- unique(  df_kastner_model_aggregate_var$roi )
roiLevels <- roiLevels[ c(9,10,11,12,13,14,7,8,1,2,3,4,5,6) ] 
df_kastner_model_aggregate_var$roi_sorted <- factor( df_kastner_model_aggregate_var$roi, levels=roiLevels ) 
par(mfrow=c(1,1))
boxplot( df_kastner_model_aggregate_var$adjR2_05 ~ df_kastner_model_aggregate_var$roi_sorted, ylim=c(0,1), frame=FALSE, las=1, col=subjCols  )
summary( lm( adjR2_05 ~ roi_sorted, data=df_kastner_model_aggregate_var) )

# all rois, kappa par
df_kastner_model_aggregate_kappa <- aggregate( kappa ~ subject * roi, mean, data= df_kastner_model_filt_clean ) 
roiLevels <- unique(  df_kastner_model_aggregate_kappa$roi )
roiLevels <- roiLevels[ c(9,10,11,12,13,14,7,8,1,2,3,4,5,6) ] 
df_kastner_model_aggregate_kappa$roi_sorted <- factor( df_kastner_model_aggregate_kappa$roi, levels=roiLevels ) 
par(mfrow=c(1,1))
boxplot( df_kastner_model_aggregate_kappa$kappa ~ df_kastner_model_aggregate_kappa$roi_sorted, ylim=c(0.2,0.8), frame=FALSE, las=1, col=subjCols  )
summary( lm( kappa ~ roi_sorted, data=df_kastner_model_aggregate_kappa) )

# all rois, AIC
df_kastner_model_aggregate_adjR2 <- aggregate( adjR2_05 ~ subject * roi, mean, data= df_kastner_model_filt_clean ) 

df_kastner_model_aggregate_adjR2_10 <- aggregate( adjR2_10 ~ subject * roi, mean, data= df_kastner_model_filt_clean ) 
df_kastner_model_aggregate_adjR2_both <- aggregate( adjR2_both ~ subject * roi, mean, data= df_kastner_model_filt_clean ) 
df_kastner_model_aggregate_adjR2_phase05 <- aggregate( adjR2_phase05 ~ subject * roi, mean, data= df_kastner_model_filt_clean ) 
df_kastner_model_aggregate_adjR2_phase10 <- aggregate( adjR2_phase10 ~ subject * roi, mean, data= df_kastner_model_filt_clean ) 

df_kastner_model_aggregate_adjR2$adjR2_10 <- df_kastner_model_aggregate_adjR2_10$adjR2_10
df_kastner_model_aggregate_adjR2$adjR2_both <- df_kastner_model_aggregate_adjR2_both$adjR2_both
df_kastner_model_aggregate_adjR2$adjR2_phase05 <- df_kastner_model_aggregate_adjR2_phase05$adjR2_phase05
df_kastner_model_aggregate_adjR2$adjR2_phase10 <- df_kastner_model_aggregate_adjR2_phase10$adjR2_phase10

roiLevels <- unique(  df_kastner_model_aggregate_adjR2$roi )
roiLevels <- roiLevels[ c(9,10,11,12,13,14,7,8,1,2,3,4,5,6) ] 
df_kastner_model_aggregate_adjR2$roi_sorted <- factor( df_kastner_model_aggregate_adjR2$roi, levels=roiLevels ) 
par(mfrow=c(1,1))
boxplot( df_kastner_model_aggregate_adjR2$adjR2_both ~ df_kastner_model_aggregate_adjR2$roi_sorted, 
         add=FALSE, boxwex=0.1, ylim=c(0, 0.5), at = 1:14-0.24, xlim=c(0,15), yaxs='i', frame=FALSE, axes=FALSE, outline=FALSE, lty=1, col='gray50'  )
boxplot( df_kastner_model_aggregate_adjR2$adjR2_05 ~ df_kastner_model_aggregate_adjR2$roi_sorted, 
         add=TRUE, frame=FALSE, las=1, boxwex=0.1, at = 1:14-0.14, xlim=c(0,15), yaxs='i', outline=FALSE, lty=1, col='lightblue' )
boxplot( df_kastner_model_aggregate_adjR2$adjR2_10 ~ df_kastner_model_aggregate_adjR2$roi_sorted, 
         add=TRUE, boxwex=0.1, at = 1:14, xlim=c(0,15), yaxs='i', frame=FALSE, axes=FALSE, outline=FALSE, lty=1, col='orange'  )
boxplot( df_kastner_model_aggregate_adjR2$adjR2_phase05 ~ df_kastner_model_aggregate_adjR2$roi_sorted, 
         add=TRUE, boxwex=0.1, at = 1:14+0.14, xlim=c(0,15), yaxs='i', frame=FALSE, axes=FALSE, outline=FALSE, lty=1, col='red'  )
boxplot( df_kastner_model_aggregate_adjR2$adjR2_phase10 ~ df_kastner_model_aggregate_adjR2$roi_sorted, 
         add=TRUE, boxwex=0.1, at = 1:14+0.24, xlim=c(0,15), yaxs='i', frame=FALSE, axes=FALSE, outline=FALSE, lty=1, col='blue'  )

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
#df_kastner_model_filt_clean$dv <- df_kastner_model_filt_clean$mu05_CCW
#df_kastner_model_filt_clean$dv <- atan2(
#  ( sin( df_kastner_model_filt_clean$mu05_CW ) + sin( df_kastner_model_filt_clean$mu05_CCW ) ) / 2, 
#  ( cos( df_kastner_model_filt_clean$mu05_CW ) + cos( df_kastner_model_filt_clean$mu05_CCW ) ) / 2
#) + pi
#df_kastner_model_filt_clean$dv <- ( df_kastner_model_filt_clean$mu05_CCW + df_kastner_model_filt_clean$mu05_CW ) / 2
df_kastner_model_filt_clean$dv <- df_kastner_model_filt_clean$phase05_CCW

for ( subjSel in unique( df_kastner_model_filt_clean$subject )  ) {
  for ( roiSel in levels( df_kastner_model_filt_clean$roi_coarse )[1:5]  ) {
    
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
