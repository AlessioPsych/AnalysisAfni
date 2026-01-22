rm( list=ls() )
source('~/abin/AFNIio.R')

mainFolder <- '/media/alessiofracasso/DATADRIVE1/Flanker'
outputFolder <- 'derivatives/resultsGlasser'
surfacesFolder <- 'derivatives/surfaceAtlases/suma_MNI152_2009_princetonAtlas'

setwd( mainFolder )
setwd( outputFolder )

dir()

#load( file='results_Glasser_MNI_Denoised.RData' )
#output_MNI_Denoised <- outputDataset
#output_MNI_Denoised$medianValue

load( file='results_Glasser_MNI_No_Denoised.RData' )
output_MNI_No_Denoised <- outputDataset
output_MNI_No_Denoised$medianValue

#load( file='results_Glasser_ORIG_Denoised.RData' )
#output_ORIG_Denoised <- outputDataset
#output_ORIG_Denoised$medianValue

#load( file='results_Glasser_ORIG_No_Denoised.RData' )
#output_ORIG_No_Denoised <- outputDataset
#output_ORIG_No_Denoised$medianValue

lhAtlas <- 'lh.Glasser_HCP_MNI.nii.gz'
rhAtlas <- 'rh.Glasser_HCP_MNI.nii.gz'

#### function definition ####

debugFlag <- 0
if (debugFlag==1) {
  sideSelected <- 'lh'
  dataSelected <- subset( output_ORIG_Denoised, side==sideSelected )
  atlasSelected <- lhAtlas
  outputFile <- 'output_Glasser_lh.nii.gz'
  fdrThreshold <- 0.01
}
generateFdrVolume <- function( sideSelected, dataSelected, atlasSelected, outputFile, fdrThreshold ) {
  unique( dataSelected$coefficient )
  str( dataSelected )
  rois <- unique( dataSelected$roi )
  testOutput <- data.frame( roi=factor(), 
                            roiIdx=double(),
                            averages=double(), 
                            stDev=double(), 
                            tStats=double(), 
                            dfs=double(),
                            storedPs=double() )
  for ( i in 1:length( rois ) ) { # i <- 1
    dataSelectedLoop <- subset( dataSelected, roi==rois[i] & coefficient=='incongruent-congruent_GLT#0_Coef' )
    tLoop <- t.test( dataSelectedLoop$medianValue )
    testOutputTemp <- data.frame( roi=dataSelectedLoop$roi[1],
                                  roiIdx=dataSelectedLoop$roiIdx[1],
                                  averages=round( median( dataSelectedLoop$medianValue ), 4 ), 
                                  stDev=round( sd( dataSelectedLoop$medianValue ), 4 ), 
                                  tStats=tLoop$statistic, 
                                  dfs=tLoop$parameter,
                                  storedPs=tLoop$p.value )
    testOutput <- rbind( testOutput, testOutputTemp )
  }
  testOutput$storePsFDR <- p.adjust( testOutput$storedPs, method=c('fdr') )
  
  significantRois <- testOutput[ testOutput$storePsFDR < fdrThreshold, ]
  
  atlasFile <- read.AFNI( atlasSelected )
  atlasVolume <- atlasFile$brk
  emptyVolume <- array(0,dim(atlasVolume))
  for ( sIdx in 1:dim( significantRois )[1] ) { #sIdx <- 1
    roiIdxTemp <- significantRois$roiIdx[sIdx]  
    whichIdxInOriginalAtlas <- which( atlasVolume==roiIdxTemp )
    emptyVolume[ whichIdxInOriginalAtlas ] <- significantRois$averages[sIdx]  
  }
  
  # delete outputfile, if it exists, output folder
  if ( file.exists( outputFile ) ) {
    instr <- sprintf( 'rm %s', outputFile )
    print( instr )
    system( instr )
  }
  # delete outputfile, if it exists, surfaces folder
  if ( file.exists( sprintf('%s/%s/%s',mainFolder,surfacesFolder,outputFile) ) ) {
    instr <- sprintf( 'rm %s', sprintf('%s/%s/%s',mainFolder,surfacesFolder,outputFile) )
    print( instr )
    system( instr )
  }
  
  write.AFNI( filename = outputFile, 
              brk=emptyVolume,
              origin=atlasFile$origin,
              defhead=atlasFile$NI_head,
              orient=atlasFile$orient )
  write.AFNI( filename = sprintf('%s/%s/%s',mainFolder,surfacesFolder,outputFile), 
              brk=emptyVolume,
              origin=atlasFile$origin,
              defhead=atlasFile$NI_head,
              orient=atlasFile$orient )
  
  return( significantRois )
  
}

#### results: output_ORIG_Denoised ####
# setwd( mainFolder )
# setwd( outputFolder )
# sideSelected <- 'lh'
# dataSelected <- subset( output_ORIG_Denoised, side==sideSelected )
# atlasSelected <- lhAtlas
# outputFileLh <- 'output_Glasser_lh_ORIG_Denoised.nii.gz'
# fdrThreshold <- 0.01
# sigRoiLh_orig_denoised <- generateFdrVolume( sideSelected, dataSelected, atlasSelected, outputFileLh, fdrThreshold )

# setwd( mainFolder )
# setwd( outputFolder )
# sideSelected <- 'rh'
# dataSelected <- subset( output_ORIG_Denoised, side==sideSelected )
# atlasSelected <- rhAtlas
# outputFileRh <- 'output_Glasser_rh_ORIG_Denoised.nii.gz'
# fdrThreshold <- 0.01
# sigRoiRh_orig_denoised <- generateFdrVolume( sideSelected, dataSelected, atlasSelected, outputFileRh, fdrThreshold )

# # clean up and combine output files
# if ( file.exists( sprintf('%s/%s/output_Glasser_ORIG_Denoised_combined.nii.gz',mainFolder,surfacesFolder) ) ) {
#   instr <- sprintf( 'rm %s', sprintf('%s/%s/output_Glasser_ORIG_Denoised_combined.nii.gz',mainFolder,surfacesFolder) )
#   print( instr )
#   system( instr )
# }
# setwd( mainFolder )
# setwd( surfacesFolder )
# instr <- sprintf( '3dcalc -a %s -b %s -expr \u0027a+b\u0027 -prefix output_Glasser_ORIG_Denoised_combined.nii.gz', outputFileLh, outputFileRh )
# print( instr )
# system( instr )
# 
# #### results: output_MNI_Denoised ####
# setwd( mainFolder )
# setwd( outputFolder )
# sideSelected <- 'lh'
# dataSelected <- subset( output_MNI_Denoised, side==sideSelected )
# atlasSelected <- lhAtlas
# outputFileLh <- 'output_Glasser_lh_MNI_Denoised.nii.gz'
# fdrThreshold <- 0.01
# sigRoiLh_MNI_denoised <- generateFdrVolume( sideSelected, dataSelected, atlasSelected, outputFileLh, fdrThreshold )
# 
# setwd( mainFolder )
# setwd( outputFolder )
# sideSelected <- 'rh'
# dataSelected <- subset( output_MNI_Denoised, side==sideSelected )
# atlasSelected <- rhAtlas
# outputFileRh <- 'output_Glasser_rh_MNI_Denoised.nii.gz'
# fdrThreshold <- 0.01
# sigRoiRh_MNI_denoised <- generateFdrVolume( sideSelected, dataSelected, atlasSelected, outputFileRh, fdrThreshold )
# 
# # clean up and combine output files
# if ( file.exists( sprintf('%s/%s/output_Glasser_MNI_Denoised_combined.nii.gz',mainFolder,surfacesFolder) ) ) {
#   instr <- sprintf( 'rm %s', sprintf('%s/%s/output_Glasser_MNI_Denoised_combined.nii.gz',mainFolder,surfacesFolder) )
#   print( instr )
#   system( instr )
# }
# setwd( mainFolder )
# setwd( surfacesFolder )
# instr <- sprintf( '3dcalc -a %s -b %s -expr \u0027a+b\u0027 -prefix output_Glasser_MNI_Denoised_combined.nii.gz', outputFileLh, outputFileRh )
# print( instr )
# system( instr )
# 
# 
# #### results: output_ORIG_No_Denoised ####
# setwd( mainFolder )
# setwd( outputFolder )
# sideSelected <- 'lh'
# dataSelected <- subset( output_ORIG_No_Denoised, side==sideSelected )
# atlasSelected <- lhAtlas
# outputFileLh <- 'output_Glasser_lh_ORIG_No_Denoised.nii.gz'
# fdrThreshold <- 0.01
# sigRoiLh_ORIG_No_denoised <- generateFdrVolume( sideSelected, dataSelected, atlasSelected, outputFileLh, fdrThreshold )
# 
# setwd( mainFolder )
# setwd( outputFolder )
# sideSelected <- 'rh'
# dataSelected <- subset(output_ORIG_No_Denoised, side==sideSelected )
# atlasSelected <- rhAtlas
# outputFileRh <- 'output_Glasser_rh_ORIG_No_Denoised.nii.gz'
# fdrThreshold <- 0.01
# sigRoiRh_ORIG_No_denoised <- generateFdrVolume( sideSelected, dataSelected, atlasSelected, outputFileRh, fdrThreshold )
# 
# # clean up and combine output files
# if ( file.exists( sprintf('%s/%s/output_Glasser_ORIG_No_Denoised_combined.nii.gz',mainFolder,surfacesFolder) ) ) {
#   instr <- sprintf( 'rm %s', sprintf('%s/%s/output_Glasser_ORIG_No_Denoised_combined.nii.gz',mainFolder,surfacesFolder) )
#   print( instr )
#   system( instr )
# }
# setwd( mainFolder )
# setwd( surfacesFolder )
# instr <- sprintf( '3dcalc -a %s -b %s -expr \u0027a+b\u0027 -prefix output_Glasser_ORIG_No_Denoised_combined.nii.gz', outputFileLh, outputFileRh )
# print( instr )
# system( instr )

#### results: output_MNI_No_Denoised ####
setwd( mainFolder )
setwd( outputFolder )
sideSelected <- 'lh'
dataSelected <- subset( output_MNI_No_Denoised, side==sideSelected )
atlasSelected <- lhAtlas
outputFileLh <- 'output_Glasser_lh_MNI_No_Denoised.nii.gz'
fdrThreshold <- 0.01
sigRoiLh_MNI_No_denoised <- generateFdrVolume( sideSelected, dataSelected, atlasSelected, outputFileLh, fdrThreshold )

setwd( mainFolder )
setwd( outputFolder )
sideSelected <- 'rh'
dataSelected <- subset(output_MNI_No_Denoised, side==sideSelected )
atlasSelected <- rhAtlas
outputFileRh <- 'output_Glasser_rh_MNI_No_Denoised.nii.gz'
fdrThreshold <- 0.01
sigRoiRh_MNI_No_denoised <- generateFdrVolume( sideSelected, dataSelected, atlasSelected, outputFileRh, fdrThreshold )

# clean up and combine output files
if ( file.exists( sprintf('%s/%s/output_Glasser_MNI_No_Denoised_combined.nii.gz',mainFolder,surfacesFolder) ) ) {
  instr <- sprintf( 'rm %s', sprintf('%s/%s/output_Glasser_MNI_No_Denoised_combined.nii.gz',mainFolder,surfacesFolder) )
  print( instr )
  system( instr )
}
setwd( mainFolder )
setwd( surfacesFolder )
instr <- sprintf( '3dcalc -a %s -b %s -expr \u0027a+b\u0027 -prefix output_Glasser_MNI_No_Denoised_combined.nii.gz', outputFileLh, outputFileRh )
print( instr )
system( instr )
