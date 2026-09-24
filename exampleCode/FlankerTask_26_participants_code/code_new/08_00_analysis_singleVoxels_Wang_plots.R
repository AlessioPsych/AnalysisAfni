rm( list=ls() )
source('~/abin/AFNIio.R')
library(lme4)
library(lmerTest)

mainFolder <- '/mnt/disk01/ds000102_FlankerTask'
outputFolder <- 'derivatives/resultsWang'
surfacesFolder <- '/mnt/disk01/surfaceAtlas_HCP_vonEconomo/surfaceAtlases/suma_MNI152_2009_princetonAtlas'

setwd( mainFolder )
setwd( outputFolder )

dir()

inputFiles <- dir( patter='*incongruent_congruent.RData' )
load( inputFiles[1] )
output_ORIG_Denoised <- outputDataset
for (i in 2:length(inputFiles) ) {
	print( sprintf('loading file %s',inputFiles[i]) )
	load( inputFiles[i] )
	dataInTemp <- outputDataset
	output_ORIG_Denoised <- rbind( output_ORIG_Denoised, dataInTemp )
}
str(output_ORIG_Denoised)

dataSelected <- output_ORIG_Denoised
dataTypeName <- 'output_ORIG_Denoised'
dataSelected$roi <- as.factor( dataSelected$roi )
setwd( mainFolder )
setwd( outputFolder  )

table( dataSelected$roi )
reorderedLevels <- levels( dataSelected$roi )[ c( 17, 19, 23, 3, 24, 25, 13, 14, 16, 18, 22, 20, 21, 10, 11, 2, 12, 4:9, 15, 1 ) ]
dataSelected$roi_recoded <- factor( dataSelected$roi, levels=reorderedLevels )
levels( dataSelected$roi_recoded )
dataSelected$macro_rois <- ifelse( dataSelected$roi_recoded==reorderedLevels[1] | 
                                     dataSelected$roi_recoded==reorderedLevels[2] |
                                     dataSelected$roi_recoded==reorderedLevels[3] |
                                     dataSelected$roi_recoded==reorderedLevels[4] |
                                     dataSelected$roi_recoded==reorderedLevels[5] |
                                     dataSelected$roi_recoded==reorderedLevels[6] |
                                     dataSelected$roi_recoded==reorderedLevels[7] |
                                     dataSelected$roi_recoded==reorderedLevels[8], 'ventral-temporal', 
                                   ifelse( dataSelected$roi_recoded==reorderedLevels[9] |
                                             dataSelected$roi_recoded==reorderedLevels[10] |
                                             dataSelected$roi_recoded==reorderedLevels[11] |
                                             dataSelected$roi_recoded==reorderedLevels[12] |
                                             dataSelected$roi_recoded==reorderedLevels[13] |
                                             dataSelected$roi_recoded==reorderedLevels[14] |
                                             dataSelected$roi_recoded==reorderedLevels[15] |
                                             dataSelected$roi_recoded==reorderedLevels[16] |
                                             dataSelected$roi_recoded==reorderedLevels[17], 'dorsal-lateral',
                                           ifelse( dataSelected$roi_recoded==reorderedLevels[18] |
                                                     dataSelected$roi_recoded==reorderedLevels[19] |
                                                     dataSelected$roi_recoded==reorderedLevels[20] |
                                                     dataSelected$roi_recoded==reorderedLevels[21] |
                                                     dataSelected$roi_recoded==reorderedLevels[22] |
                                                     dataSelected$roi_recoded==reorderedLevels[23] |
                                                     dataSelected$roi_recoded==reorderedLevels[24] |
                                                     dataSelected$roi_recoded==reorderedLevels[25], 'parietal-frontal', 'none' ) ) )
table( dataSelected$macro_rois, dataSelected$roi_recoded )
dataSelected_agg <- aggregate( voxelValue ~ participantID * roi_recoded * coefficient * macro_rois, dataSelected, FUN=median )
dataSelected_agg$medianValue <- dataSelected_agg$voxelValue

x11( width=8.5, height=7 )
par(mfrow=c(2,2))

#### plot ventral - temporal ####

dataPlot <- subset( dataSelected_agg, macro_rois=='ventral-temporal' ) 
dataPlot_cong <- subset( dataPlot, coefficient=='congruent#0_Coef' )
dataPlot_cong$roi_recoded <- droplevels( dataPlot_cong$roi_recoded )
dataPlot_incong <- subset( dataPlot, coefficient=='incongruent#0_Coef' )
dataPlot_incong$roi_recoded <- droplevels( dataPlot_incong$roi_recoded )

boxplot( dataPlot_cong$medianValue ~ dataPlot_cong$roi_recoded, 
         outline=FALSE, frame=FALSE, las=1, col='lightblue',
         ylim=c(-0.4,1), xlim=c(0,8.5), at=seq(0.4,7.5,1), 
         axes=FALSE, boxwex=0.3, xlab='', ylab='' )
boxplot( dataPlot_incong$medianValue ~ dataPlot_cong$roi_recoded, 
         outline=FALSE, frame=FALSE, las=1, col='orange',
         ylim=c(-0.4,1), xlim=c(0,8.5), at=seq(0.7,7.7,1), 
         axes=FALSE, boxwex=0.3, add=TRUE )
axis(1, seq(0.5,7.6,1), levels( dataPlot_cong$roi_recoded ), las=2, cex.axis=1.4 )
axis(2, c(-0.4, 0, 0.5, 1), c(-0.4, 0, 0.5, 1), las=1, cex.axis=1.4 )
abline( h=0, lwd=2, lty=2 )
roiLevels <- levels( dataPlot_cong$roi_recoded )
xLabelPosition <- seq(0.5,7.6,1)
for (k in 1:length( roiLevels ) ) { #k <- 1
  sampleCong <- dataPlot_cong$medianValue[ dataPlot_cong$roi_recoded==roiLevels[k] ]  
  sampleInCong <- dataPlot_incong$medianValue[ dataPlot_incong$roi_recoded==roiLevels[k] ]  
  storeT_avgSignal <- t.test( sampleCong + sampleInCong ) # test average signal
  storeT <- t.test( sampleCong, sampleInCong, paired=TRUE ) # test difference between conditions
  corrPavgsignal <- storeT_avgSignal$p.value*25
  corrPconditions <- storeT$p.value*25
  
  print( sprintf( '%s: %1.4f avg signal', roiLevels[k], corrPavgsignal ) )
  print( sprintf( '%s: %1.4f diff conditions', roiLevels[k], corrPconditions ) )
  if ( corrPavgsignal<0.05 ) {
    points( xLabelPosition[k], 0.9, type='p', pch=0, cex=1.5, las=2 )
    if ( corrPconditions<0.05 ) {
      points( xLabelPosition[k], 1, type='p', pch=2, cex=1.5, las=2 )
    }
  }
}

# interaction / main effect approach
depVariable <- c( dataPlot_cong$medianValue, dataPlot_incong$medianValue )
indVar_ROI <- as.factor( c( dataPlot_cong$roi_recoded, dataPlot_incong$roi_recoded ) )
indVar_congVar <- as.factor( c( rep('cong', length( dataPlot_cong$roi_recoded ) ), 
                                rep('incong', length( dataPlot_incong$roi_recoded ) ) ) )
partVar <- as.factor( c( dataPlot_cong$participantID, dataPlot_incong$participantID ) )
table( partVar, indVar_ROI, indVar_congVar )
mod <- lmer( depVariable ~ indVar_ROI * indVar_congVar + (1|partVar) )
summary( mod )
anova( mod )

modAov <- aov( depVariable ~ indVar_ROI * indVar_congVar + Error(partVar) )
summary( modAov )

#### plot dorsal - lateral ####

dataPlot <- subset( dataSelected_agg, macro_rois=='dorsal-lateral' ) 
dataPlot_cong <- subset( dataPlot, coefficient=='congruent#0_Coef' )
dataPlot_cong$roi_recoded <- droplevels( dataPlot_cong$roi_recoded )
dataPlot_incong <- subset( dataPlot, coefficient=='incongruent#0_Coef' )
dataPlot_incong$roi_recoded <- droplevels( dataPlot_incong$roi_recoded )

boxplot( dataPlot_cong$medianValue ~ dataPlot_cong$roi_recoded, 
         outline=FALSE, frame=FALSE, las=1, col='lightblue',
         ylim=c(-0.4,1), xlim=c(0,9.5), at=seq(0.4,8.5,1), 
         axes=FALSE, boxwex=0.3, xlab='', ylab='' )
boxplot( dataPlot_incong$medianValue ~ dataPlot_cong$roi_recoded, 
         outline=FALSE, frame=FALSE, las=1, col='orange',
         ylim=c(-0.4,1), xlim=c(0,9.5), at=seq(0.7,8.7,1), 
         axes=FALSE, boxwex=0.3, add=TRUE )
axis(1, seq(0.5,8.6,1), levels( dataPlot_cong$roi_recoded ), las=2, cex.axis=1.4 )
axis(2, c(-0.4, 0, 0.5, 1), c(-0.4, 0, 0.5, 1), las=1, cex.axis=1.4 )
abline( h=0, lwd=2, lty=2 )
roiLevels <- levels( dataPlot_cong$roi_recoded )
xLabelPosition <- seq(0.5,8.6,1)
for (k in 1:length( roiLevels ) ) { #k <- 1
  sampleCong <- dataPlot_cong$medianValue[ dataPlot_cong$roi_recoded==roiLevels[k] ]  
  sampleInCong <- dataPlot_incong$medianValue[ dataPlot_incong$roi_recoded==roiLevels[k] ]  
  storeT_avgSignal <- t.test( sampleCong + sampleInCong ) # test average signal
  storeT <- t.test( sampleCong, sampleInCong, paired=TRUE ) # test difference between conditions
  corrPavgsignal <- storeT_avgSignal$p.value*25
  corrPconditions <- storeT$p.value*25
  
  print( sprintf( '%s: %1.4f avg signal', roiLevels[k], corrPavgsignal ) )
  print( sprintf( '%s: %1.4f diff conditions', roiLevels[k], corrPconditions ) )
  if ( corrPavgsignal<0.05 ) {
    points( xLabelPosition[k], 0.9, type='p', pch=0, cex=1.5, las=2 )
    if ( corrPconditions<0.05 ) {
      points( xLabelPosition[k], 1, type='p', pch=2, cex=1.5, las=2 )
    }
  }
}

# interaction / main effect approach
depVariable <- c( dataPlot_cong$medianValue, dataPlot_incong$medianValue )
indVar_ROI <- as.factor( c( dataPlot_cong$roi_recoded, dataPlot_incong$roi_recoded ) )
indVar_congVar <- as.factor( c( rep('cong', length( dataPlot_cong$roi_recoded ) ), 
                                rep('incong', length( dataPlot_incong$roi_recoded ) ) ) )
partVar <- as.factor( c( dataPlot_cong$participantID, dataPlot_incong$participantID ) )
table( partVar, indVar_ROI, indVar_congVar )
mod <- lmer( depVariable ~ indVar_ROI * indVar_congVar + (1|partVar) )
summary( mod )
anova( mod )

modAov <- aov( depVariable ~ indVar_ROI * indVar_congVar + Error(partVar) )
summary( modAov )


#### plot parietal - frontal ####

dataPlot <- subset( dataSelected_agg, macro_rois=='parietal-frontal' ) 
dataPlot_cong <- subset( dataPlot, coefficient=='congruent#0_Coef' )
dataPlot_cong$roi_recoded <- droplevels( dataPlot_cong$roi_recoded )
dataPlot_incong <- subset( dataPlot, coefficient=='incongruent#0_Coef' )
dataPlot_incong$roi_recoded <- droplevels( dataPlot_incong$roi_recoded )

boxplot( dataPlot_cong$voxelValue ~ dataPlot_cong$roi_recoded, 
         outline=FALSE, frame=FALSE, las=1, col='lightblue',
         ylim=c(-0.4,1), xlim=c(0,8.5), at=seq(0.4,7.5,1), 
         axes=FALSE, boxwex=0.3, xlab='', ylab='' )
boxplot( dataPlot_incong$voxelValue ~ dataPlot_cong$roi_recoded, 
         outline=FALSE, frame=FALSE, las=1, col='orange',
         ylim=c(-0.4,1), xlim=c(0,8.5), at=seq(0.7,7.7,1), 
         axes=FALSE, boxwex=0.3, add=TRUE )
axis(1, seq(0.5,7.6,1), levels( dataPlot_cong$roi_recoded ), las=2, cex.axis=1.4 )
axis(2, c(-0.4, 0, 0.5, 1), c(-0.4, 0, 0.5, 1), las=1, cex.axis=1.4 )
abline( h=0, lwd=2, lty=2 )
roiLevels <- levels( dataPlot_cong$roi_recoded )
xLabelPosition <- seq(0.5,7.6,1)
for (k in 1:length( roiLevels ) ) { #k <- 1
  sampleCong <- dataPlot_cong$medianValue[ dataPlot_cong$roi_recoded==roiLevels[k] ]  
  sampleInCong <- dataPlot_incong$medianValue[ dataPlot_incong$roi_recoded==roiLevels[k] ]  
  storeT_avgSignal <- t.test( sampleCong + sampleInCong ) # test average signal
  storeT <- t.test( sampleCong, sampleInCong, paired=TRUE ) # test difference between conditions
  corrPavgsignal <- storeT_avgSignal$p.value*25
  corrPconditions <- storeT$p.value*25
  
  print( sprintf( '%s: %1.4f avg signal', roiLevels[k], corrPavgsignal ) )
  print( sprintf( '%s: %1.4f diff conditions', roiLevels[k], corrPconditions ) )
  if ( corrPavgsignal<0.05 ) {
    points( xLabelPosition[k], 0.9, type='p', pch=0, cex=1.5, las=2 )
    if ( corrPconditions<0.05 ) {
      points( xLabelPosition[k], 1, type='p', pch=2, cex=1.5, las=2 )
    }
  }
}

#### trend along the intra-parietal sulcus
dataPlot_cong$cong_incong_diff <- dataPlot_incong$medianValue - dataPlot_cong$medianValue
dataPlotTest_IPS <- subset( dataPlot_cong, dataPlot_cong$roi_recoded=='IPS0' | 
                              dataPlot_cong$roi_recoded=='IPS1' | 
                              dataPlot_cong$roi_recoded=='IPS2' |
                              dataPlot_cong$roi_recoded=='IPS3' |
                              dataPlot_cong$roi_recoded=='IPS4' |
                              dataPlot_cong$roi_recoded=='IPS5' )
dataPlotTest_IPS$cong_incong_diff
dataPlotTest_IPS$roi_recoded <- droplevels( dataPlotTest_IPS$roi_recoded )
levels( dataPlotTest_IPS$roi_recoded )
table( dataPlotTest_IPS$roi_recoded )
dataPlotTest_IPS$roi_recoded_numeric <- as.numeric( dataPlotTest_IPS$roi_recoded )
boxplot( dataPlotTest_IPS$cong_incong_diff ~ dataPlotTest_IPS$roi_recoded, outline=FALSE,
         frame=FALSE, las=1, xlab='', ylab='% BOLD signal difference', ylim=c(-0.2,0.7),
         cex.lab=1.4, cex.axis=1.4, col='grey80', boxwex=0.5, axes=FALSE )
axis( 1, c(1,2,3,4,5,6), c('IPS0','IPS1','IPS2','IPS3','IPS4','IPS5'), las=2, cex.axis=1.4 )
axis( 2, c(0,0.25,0.5),c(0,0.25,0.5), las=1, cex.axis=1.4 )
mod <- lm( dataPlotTest_IPS$cong_incong_diff ~ dataPlotTest_IPS$roi_recoded_numeric )
abline( mod, col='grey50', lwd=4, lty=2 )
summary( mod )

#### trend along the intra-parietal sulcus, interaction approach
library( lme4 )
library( lmerTest )
dataPlot_cong$incong_medianValue <- dataPlot_incong$medianValue
dataPlotTest_IPS <- subset( dataPlot_cong, dataPlot_cong$roi_recoded=='IPS0' | 
                              dataPlot_cong$roi_recoded=='IPS1' | 
                              dataPlot_cong$roi_recoded=='IPS2' |
                              dataPlot_cong$roi_recoded=='IPS3' |
                              dataPlot_cong$roi_recoded=='IPS4' |
                              dataPlot_cong$roi_recoded=='IPS5' )
dataPlotTest_IPS$incong_medianValue
dataPlotTest_IPS$roi_recoded <- droplevels( dataPlotTest_IPS$roi_recoded )
levels( dataPlotTest_IPS$roi_recoded )
table( dataPlotTest_IPS$roi_recoded )
dataPlotTest_IPS$roi_recoded_numeric <- as.numeric( dataPlotTest_IPS$roi_recoded )

depVariable <- c( dataPlotTest_IPS$medianValue, dataPlotTest_IPS$incong_medianValue )
indVar_ROI <- as.factor( c( dataPlotTest_IPS$roi_recoded, dataPlotTest_IPS$roi_recoded ) )
indVar_congVar <- as.factor( c( rep('cong', length( dataPlotTest_IPS$roi_recoded ) ), 
                                rep('incong', length( dataPlotTest_IPS$roi_recoded ) ) ) )
partVar <- as.factor( c( dataPlotTest_IPS$participantID, dataPlotTest_IPS$participantID ) )
table( partVar, indVar_ROI, indVar_congVar )
mod <- lmer( depVariable ~ indVar_ROI * indVar_congVar + (1|partVar) )
modAov <- aov( depVariable ~ indVar_ROI * indVar_congVar + Error(partVar) )
summary( mod )
anova( mod )
summary( modAov )

dev.copy2pdf( file=sprintf('%s_resultsBOLD.pdf',dataTypeName), width=8.5, height=7 )
dev.off()

x11( width=8.5, height=7 )
par(mfrow=c(2,2))

#### plot parietal - frontal, tsnr ####

dataPlot <- subset( dataSelected_agg, macro_rois=='parietal-frontal' ) 
dataPlot_cong <- subset( dataPlot, coefficient=='tsnr' )
dataPlot_cong$roi_recoded <- droplevels( dataPlot_cong$roi_recoded )
dataPlot_incong <- subset( dataPlot, coefficient=='tsnr' )
dataPlot_incong$roi_recoded <- droplevels( dataPlot_incong$roi_recoded )

boxplot( dataPlot_cong$medianValue ~ dataPlot_cong$roi_recoded, 
         outline=FALSE, frame=FALSE, las=1, col='lightblue',
         ylim=c(100,400), xlim=c(0,8.5), at=seq(0.4,7.5,1), 
         axes=FALSE, boxwex=0.8, xlab='', ylab='' )
axis(1, seq(0.5,7.6,1), levels( dataPlot_cong$roi_recoded ), las=2, cex.axis=1.4 )
axis(2, c(100, 200, 300, 400), c(100, 200, 300, 400), las=1, cex.axis=1.4 )

#### plot ventral - temporal ####

dataPlot <- subset( dataSelected_agg, macro_rois=='ventral-temporal' ) 
dataPlot_cong <- subset( dataPlot, coefficient=='tsnr' )
dataPlot_cong$roi_recoded <- droplevels( dataPlot_cong$roi_recoded )
dataPlot_incong <- subset( dataPlot, coefficient=='tsnr' )
dataPlot_incong$roi_recoded <- droplevels( dataPlot_incong$roi_recoded )

boxplot( dataPlot_cong$medianValue ~ dataPlot_cong$roi_recoded, 
         outline=FALSE, frame=FALSE, las=1, col='lightblue',
         ylim=c(100,400), xlim=c(0,8.5), at=seq(0.4,7.5,1), 
         axes=FALSE, boxwex=0.8, xlab='', ylab='' )
axis(1, seq(0.5,7.6,1), levels( dataPlot_cong$roi_recoded ), las=2, cex.axis=1.4 )
axis(2, c(100, 200, 300, 400), c(100, 200, 300, 400), las=1, cex.axis=1.4 )

#### plot ventral - temporal ####

dataPlot <- subset( dataSelected_agg, macro_rois=='dorsal-lateral' ) 
dataPlot_cong <- subset( dataPlot, coefficient=='tsnr' )
dataPlot_cong$roi_recoded <- droplevels( dataPlot_cong$roi_recoded )
dataPlot_incong <- subset( dataPlot, coefficient=='tsnr' )
dataPlot_incong$roi_recoded <- droplevels( dataPlot_incong$roi_recoded )

boxplot( dataPlot_cong$medianValue ~ dataPlot_cong$roi_recoded, 
         outline=FALSE, frame=FALSE, las=1, col='lightblue',
         ylim=c(100,400), xlim=c(0,9.5), at=seq(0.4,8.5,1), 
         axes=FALSE, boxwex=0.8, xlab='', ylab='' )
axis(1, seq(0.5,8.6,1), levels( dataPlot_cong$roi_recoded ), las=2, cex.axis=1.4 )
axis(2, c(100, 200, 300, 400), c(100, 200, 300, 400), las=1, cex.axis=1.4 )

dev.copy2pdf( file=sprintf('%s_resultsTSNR.pdf',dataTypeName), width=8.5, height=7 )
dev.off()



















