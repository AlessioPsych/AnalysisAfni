rm(list=ls()); gc();

# toolbox and functions
generalPurposeDir <- Sys.getenv( x='AFNI_TOOLBOXDIRGENERALPURPOSE' )
afniInstallDir <- Sys.getenv( x='AFNI_INSTALLDIR' ) 
source( sprintf('%s/scaleData.R', generalPurposeDir) )
source( sprintf('%s/AFNIio.R', afniInstallDir ) )
library( pracma )
library( abind )
library( circular )

# folders 
mainFolder <- '/media/alessiofracasso/DATADRIVE1/KastnerData/staging_area_Kastner'
anatomyDir <- sprintf( '%s/anatomies_KastnerClassic_Freesurfer', mainFolder )
modelOutputFolder <- sprintf('%s/modelsOutput/kastnerClassic_contrastBased', mainFolder )
modelOutputFolder_kastnerModel <- sprintf('%s/modelsOutput/kastnerClassic', mainFolder )
surfaceFolder <- sprintf('%s/anatomies_KastnerClassic_Freesurfer/surfaceAtlases/suma_MNI152_2009', mainFolder)
roiSurfaceFolderMNI <- sprintf('%s/probatlas_v4/ProbAtlas_v4/subj_surf_all', surfaceFolder)
resultsFolder <- sprintf('%s/results', mainFolder ) 
figureFolder <- sprintf('%s/figures_KastnerClassic_hrf_modulation',mainFolder)

setwd( modelOutputFolder )
print( sprintf('current folder: %s', getwd() ) )
modelParticipants <- dir( modelOutputFolder )
print( sprintf('model participants:' ) )
modelParticipants

setwd( modelOutputFolder_kastnerModel )
print( sprintf('current folder: %s', getwd() ) )
modelParticipants_kastnerModel <- dir( modelOutputFolder_kastnerModel )
print( sprintf('model participants kastner classic HRF:' ) )
modelParticipants_kastnerModel

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
anatomyFolders_matching <- rep(0,length(modelParticipants))
for ( nPart in 1:length( modelParticipants ) ) {
  partId <- strsplit( modelParticipants[ nPart ], '_' )[[1]][1]
  anatIdx <- which( anatomyFolders_id == partId )
  print( anatIdx )
  anatomyFolders_matching[nPart] <- anatIdx
}

# check the alignment
selectedAnatFolders <- anatomyFolders[ anatomyFolders_matching ]
data.frame( modelParticipants, modelParticipants_kastnerModel, selectedAnatFolders )

runCodeFlag <- 1
Sys.setenv(OMP_NUM_THREADS='8')
Sys.getenv('OMP_NUM_THREADS')

##########################
### select participant ###
##########################
### select roi, model, hemi ###
###############################
nPart <- 7
selectedRoiIdx <- 21
selectedHemiRoisFlag <- 'left'
selectedModelFlag <- 'cw'
##########################

# move to the right anatomy folder, where modelling results are stored
setwd( anatomyDir )
setwd( selectedAnatFolders[ nPart ] )
setwd( 'FreeSeg_result/SUMA' )
setwd( 'roisWang2015_KastnerClassic' )
print('########')
print( getwd() )
print('########')

# get roi data:
if ( selectedHemiRoisFlag=='left' ) {
  roi <- read.AFNI( filename='WangRoi2015.lh.nii.gz' )
  selectedHemiRois <- roi$brk[,,,1]
}
if ( selectedHemiRoisFlag=='right' ) {
  roi <- read.AFNI( filename='WangRoi2015.rh.nii.gz' )
  selectedHemiRois <- roi$brk[,,,1]
}
rm( list=c('roi') ); gc(); gc() #clean up

#############################################
### load model kastner prf contrast based ###
#############################################

setwd( modelOutputFolder )
print( sprintf('current folder: %s', getwd() ) )
setwd( modelParticipants[ nPart ] )
print('########')
print( getwd() )
print('########')

# get model data:
if ( selectedModelFlag=='cw' ) {
  modelContrast <- read.AFNI( filename=sprintf('%s_cw_R2_pRFBasedPrediction.nii.gz', modelParticipants[ nPart ] ) )
  selectedModel <- modelContrast$brk
  modelContrastTs <- read.AFNI( filename=sprintf('%s_cw_detrendedTs_pRFBasedPrediction.nii.gz', modelParticipants[ nPart ] ) )
  selectedTsDetrended <- modelContrastTs$brk
  modelContrastPredictedTs <- read.AFNI( filename=sprintf('%s_cw_predictedTs_pRFBasedPrediction.nii.gz', modelParticipants[ nPart ] ) )
  selectedTsPredicted <- modelContrastPredictedTs$brk
}
if ( selectedModelFlag=='ccw' ) {
  modelContrast <- read.AFNI( filename=sprintf('%s_ccw_R2_pRFBasedPrediction.nii.gz', modelParticipants[ nPart ] ) )
  selectedModel <- modelContrast$brk
  modelContrastTs <- read.AFNI( filename=sprintf('%s_ccw_detrendedTs_pRFBasedPrediction.nii.gz', modelParticipants[ nPart ] ) )
  selectedTsDetrended <- modelContrastTs$brk
  modelContrastPredictedTs <- read.AFNI( filename=sprintf('%s_ccw_predictedTs_pRFBasedPrediction.nii.gz', modelParticipants[ nPart ] ) )
  selectedTsPredicted <- modelContrastPredictedTs$brk
}
rm( list=c('modelContrast','modelContrastTs','modelContrastPredictedTs') ); gc(); gc() #clean up


#############################################
### load model kastner classic, 5 cycles ####
#############################################

setwd( modelOutputFolder_kastnerModel )
print( sprintf('current folder: %s', getwd() ) )
setwd( modelParticipants[ nPart ] )
print('########')
print( getwd() )
print('########')

# get model data:
if ( selectedModelFlag=='cw' ) {
  modelContrast <- read.AFNI( filename=sprintf('%s_CW_2components_05_params.nii.gz', modelParticipants_kastnerModel[ nPart ] ) )
  selectedModel_kastner <- modelContrast$brk
  modelContrastTs <- read.AFNI( filename=sprintf('%s_CW_2components_05_DetrendedTs.nii.gz', modelParticipants[ nPart ] ) )
  selectedTsDetrended_kastner <- modelContrastTs$brk
  modelContrastPredictedTs <- read.AFNI( filename=sprintf('%s_CW_2components_05_PredixtedTs.nii.gz', modelParticipants[ nPart ] ) )
  selectedTsPredicted_kastner <- modelContrastPredictedTs$brk
}
if ( selectedModelFlag=='ccw' ) {
  modelContrast <- read.AFNI( filename=sprintf('%s_CCW_2components_05_params.nii.gz', modelParticipants_kastnerModel[ nPart ] ) )
  selectedModel_kastner <- modelContrast$brk
  modelContrastTs <- read.AFNI( filename=sprintf('%s_CCW_2components_05_DetrendedTs.nii.gz', modelParticipants[ nPart ] ) )
  selectedTsDetrended_kastner <- modelContrastTs$brk
  modelContrastPredictedTs <- read.AFNI( filename=sprintf('%s_CCW_2components_05_PredixtedTs.nii.gz', modelParticipants[ nPart ] ) )
  selectedTsPredicted_kastner <- modelContrastPredictedTs$brk
}
rm( list=c('modelContrast','modelContrastTs','modelContrastPredictedTs') ); gc(); gc() #clean up

modelMatArray <- array( selectedModel, c( prod( dim( selectedModel )[1:3] ), dim( selectedModel )[4] ) )
detrendedTsMatArray <- array( selectedTsDetrended, c( prod( dim( selectedTsDetrended )[1:3] ), dim( selectedTsDetrended )[4] ) )
predictedTsMatArray <- array( selectedTsPredicted, c( prod( dim( selectedTsPredicted )[1:3] ), dim( selectedTsPredicted )[4] ) )

modelMatArray_kastner <- array( selectedModel_kastner, c( prod( dim( selectedModel_kastner )[1:3] ), dim( selectedModel_kastner )[4] ) )
detrendedTsMatArray_kastner <- array( selectedTsDetrended_kastner, c( prod( dim( selectedTsDetrended_kastner )[1:3] ), dim( selectedTsDetrended_kastner )[4] ) )
predictedTsMatArray_kastner <- array( selectedTsPredicted_kastner, c( prod( dim( selectedTsPredicted_kastner )[1:3] ), dim( selectedTsPredicted_kastner )[4] ) )

roiArray <- array( selectedHemiRois, c( prod( dim( selectedHemiRois )[1:3] ), 1 ) )

r2Model <- modelMatArray[,1]
r2Model_kastner <- modelMatArray_kastner[,14]

r2ModelRoi <- r2Model[ roiArray==selectedRoiIdx ]
detrendedTsMatArrayRoi <- detrendedTsMatArray[ roiArray==selectedRoiIdx, ]
predictedTsMatArrayRoi <- predictedTsMatArray[ roiArray==selectedRoiIdx, ]
r2ModelRoi_kastner <- r2Model_kastner[ roiArray==selectedRoiIdx ]
detrendedTsMatArrayRoi_kastner <- detrendedTsMatArray_kastner[ roiArray==selectedRoiIdx, ]
predictedTsMatArrayRoi_kastner <- predictedTsMatArray_kastner[ roiArray==selectedRoiIdx, ]

# model roi based on cotrast based responses OR kastner model
#r2ModelSorted <- sort( r2ModelRoi_kastner, decreasing=TRUE, index.return=TRUE )
for ( nPlot in 1:9 ) {
  r2ModelSorted <- sort( r2ModelRoi, decreasing=TRUE, index.return=TRUE )
  r2ModelSorted$ix
  selectedIdx <- nPlot
  selectedVoxelIdx <- r2ModelSorted$ix[ selectedIdx ]
  
  r2Voxel <- r2ModelRoi[ selectedVoxelIdx ]
  detrendVoxel <- detrendedTsMatArrayRoi[ selectedVoxelIdx, ]
  predictedVoxel <- predictedTsMatArrayRoi[ selectedVoxelIdx, ]
  r2Voxel_kastner <- r2ModelRoi_kastner[ selectedVoxelIdx ]
  detrendVoxel_kastner <- detrendedTsMatArrayRoi_kastner[ selectedVoxelIdx, ]
  predictedVoxel_kastner <- predictedTsMatArrayRoi_kastner[ selectedVoxelIdx, ]
  timeVar <- seq( 0, length( detrendVoxel ), length.out=length( detrendVoxel ) ) * 2
  
  
  
  # for reference, roi idxs, wang atlas
  roiNames <- c('V1v','V1d','V2v','V2d','V3v','V3d','hV4','VO1','VO2','PHC1','PHC2','MST','hMT',	   
                'LO2','LO1','V3b','V3a','IPS0','IPS1','IPS2','IPS3','IPS4','IPS5','SPL1','FEF')	
  filenameFigure <- sprintf('%s/single_time_series_SortedContrast/subj_%s_area_%s_roi_%s_orientation_%s_plot_%d.png', 
                            figureFolder,
                            selectedAnatFolders[ nPart ],
                            roiNames[selectedRoiIdx],
                            selectedHemiRoisFlag,
                            selectedModelFlag,
                            nPlot )
  widthPar <- 11
  heightPar <- 4
  x11( width=widthPar, height=heightPar )
  
  png( filename=filenameFigure, width=widthPar, height=heightPar, units='in', 
       pointsize=12, res=150, type=c('cairo') )
  
  par(mfrow=c(1,2))
  par( mar=c(5,5,5,5) )
  par(mgp=c(2.5,1,0))
  #detrendVoxelPlot <- detrendVoxel
  #predictedVoxelPlot <- predictedVoxel
  detrendVoxelPlot <- detrendVoxel_kastner
  modLmContrast <- lm( detrendVoxelPlot ~ predictedVoxel )
  predictedVoxelPlot <- predict( modLmContrast )
  r2Voxel <- summary( modLmContrast )$adj.r.squared
  plot( detrendVoxelPlot ~ timeVar, bty='n', type='b', pch=16,
        cex=1.25, col='grey70', xlab='time (s)', ylab='BOLD signal (A.U.)',
        las=1, cex.lab=1.5, cex.axis=1.5, lwd=1, lty=2 )
  lines( timeVar, predictedVoxelPlot, col='darkorange', lwd=2, lty=1 )
  lines( timeVar, predictedVoxelPlot, col='darkorange', lwd=2, lty=1 )
  text( 250, min( detrendVoxelPlot ), 
        sprintf( 'adj.R2=%1.2f', round( r2Voxel, 4 ) ), cex=1.25 )  
  summary(lm(detrendVoxelPlot~predictedVoxelPlot))
  
  par( mar=c(5,5,5,5) )
  par(mgp=c(2.5,1,0))
  detrendVoxelPlot <- detrendVoxel_kastner
  predictedVoxelPlot <- predictedVoxel_kastner
  plot( detrendVoxelPlot ~ timeVar, bty='n', type='b', pch=16,
        cex=1.25, col='grey70', xlab='time (s)', ylab='BOLD signal (A.U.)',
        las=1, cex.lab=1.5, cex.axis=1.5, lwd=1, lty=2 )
  lines( timeVar, predictedVoxelPlot, col='darkorange', lwd=2, lty=1 )
  lines( timeVar, predictedVoxelPlot, col='darkorange', lwd=2, lty=1 )
  text( 250, min( detrendVoxelPlot ), 
        sprintf( 'adj.R2=%1.2f', round( r2Voxel_kastner, 2 ) ), cex=1.25 )  
  summary(lm(detrendVoxelPlot~predictedVoxelPlot))
  
  
  #dev.copy(device=x11)
  #savePlot( filename=filenameFigure, type = c('png'))
  #dev.copy2pdf( file=filenameFigure, width=widthPar, height=heightPar )
  dev.off()
  graphics.off()
}

for ( nPlot in 1:9 ) {
  r2ModelSorted <- sort( r2ModelRoi_kastner, decreasing=TRUE, index.return=TRUE )
  r2ModelSorted$ix
  selectedIdx <- nPlot
  selectedVoxelIdx <- r2ModelSorted$ix[ selectedIdx ]
  
  r2Voxel <- r2ModelRoi[ selectedVoxelIdx ]
  detrendVoxel <- detrendedTsMatArrayRoi[ selectedVoxelIdx, ]
  predictedVoxel <- predictedTsMatArrayRoi[ selectedVoxelIdx, ]
  r2Voxel_kastner <- r2ModelRoi_kastner[ selectedVoxelIdx ]
  detrendVoxel_kastner <- detrendedTsMatArrayRoi_kastner[ selectedVoxelIdx, ]
  predictedVoxel_kastner <- predictedTsMatArrayRoi_kastner[ selectedVoxelIdx, ]
  timeVar <- seq( 0, length( detrendVoxel ), length.out=length( detrendVoxel ) ) * 2
  
  # for reference, roi idxs, wang atlas
  roiNames <- c('V1v','V1d','V2v','V2d','V3v','V3d','hV4','VO1','VO2','PHC1','PHC2','MST','hMT',	   
                'LO2','LO1','V3b','V3a','IPS0','IPS1','IPS2','IPS3','IPS4','IPS5','SPL1','FEF')	
  filenameFigure <- sprintf('%s/single_time_series_SortedKastner/subj_%s_area_%s_roi_%s_orientation_%s_plot_%d.png', 
                            figureFolder,
                            selectedAnatFolders[ nPart ],
                            roiNames[selectedRoiIdx],
                            selectedHemiRoisFlag,
                            selectedModelFlag,
                            nPlot )
  widthPar <- 11
  heightPar <- 4
  x11( width=widthPar, height=heightPar )
  
  png( filename=filenameFigure, width=widthPar, height=heightPar, units='in', 
       pointsize=12, res=150, type=c('cairo') )
  
  par(mfrow=c(1,2))
  par( mar=c(5,5,5,5) )
  par(mgp=c(2.5,1,0))
  #detrendVoxelPlot <- detrendVoxel
  #predictedVoxelPlot <- predictedVoxel
  detrendVoxelPlot <- detrendVoxel_kastner
  modLmContrast <- lm( detrendVoxelPlot ~ predictedVoxel )
  predictedVoxelPlot <- predict( modLmContrast )
  r2Voxel <- summary( modLmContrast )$adj.r.squared
  plot( detrendVoxelPlot ~ timeVar, bty='n', type='b', pch=16,
        cex=1.25, col='grey70', xlab='time (s)', ylab='BOLD signal (A.U.)',
        las=1, cex.lab=1.5, cex.axis=1.5, lwd=1, lty=2 )
  lines( timeVar, predictedVoxelPlot, col='darkorange', lwd=2, lty=1 )
  lines( timeVar, predictedVoxelPlot, col='darkorange', lwd=2, lty=1 )
  text( 250, min( detrendVoxelPlot ), 
        sprintf( 'adj.R2=%1.2f', round( r2Voxel, 4 ) ), cex=1.25 )  
  summary(lm(detrendVoxelPlot~predictedVoxelPlot))
  
  par( mar=c(5,5,5,5) )
  par(mgp=c(2.5,1,0))
  detrendVoxelPlot <- detrendVoxel_kastner
  predictedVoxelPlot <- predictedVoxel_kastner
  plot( detrendVoxelPlot ~ timeVar, bty='n', type='b', pch=16,
        cex=1.25, col='grey70', xlab='time (s)', ylab='BOLD signal (A.U.)',
        las=1, cex.lab=1.5, cex.axis=1.5, lwd=1, lty=2 )
  lines( timeVar, predictedVoxelPlot, col='darkorange', lwd=2, lty=1 )
  lines( timeVar, predictedVoxelPlot, col='darkorange', lwd=2, lty=1 )
  text( 250, min( detrendVoxelPlot ), 
        sprintf( 'adj.R2=%1.2f', round( r2Voxel_kastner, 2 ) ), cex=1.25 )  
  summary(lm(detrendVoxelPlot~predictedVoxelPlot))
  
  
  #dev.copy(device=x11)
  #savePlot( filename=filenameFigure, type = c('png'))
  #dev.copy2pdf( file=filenameFigure, width=widthPar, height=heightPar )
  dev.off()
  graphics.off()
}

#nPart <- 3
#selectedRoiIdx <- 1
#selectedHemiRoisFlag <- 'left'
#selectedModelFlag <- 'cw'










