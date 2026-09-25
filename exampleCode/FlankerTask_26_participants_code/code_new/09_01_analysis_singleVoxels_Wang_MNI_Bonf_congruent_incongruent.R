rm( list=ls() )
source('~/abin/AFNIio.R')
library(lme4)
library(lmerTest)

# to run ' Rscript filename '


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
dataSelected_agg <- aggregate( voxelValue ~ participantID * roi * roiIdx * side * coefficient, dataSelected, FUN=median )
dataSelected_agg$medianValue <- dataSelected_agg$voxelValue


lhAtlas <- 'lh.Wang_2015.nii.gz'
rhAtlas <- 'rh.Wang_2015.nii.gz'

#### function definition ####

debugFlag <- 0
if (debugFlag==1) {
  sideSelected <- 'lh'
  dataSelected <- subset( dataSelected_agg, side==sideSelected )
  atlasSelected <- lhAtlas
  outputFile <- 'output_Wang_lh.nii.gz'
  fdrThreshold <- 0.05
}

generateFdrVolume <- function( sideSelected, dataSelected, atlasSelected, outputFile, fdrThreshold, surfacesFolder ) {
  unique( dataSelected$coefficient )
  str( dataSelected )
  rois <- unique( dataSelected$roi )
  testOutput <- data.frame( roi=factor(),
  			    roiIdx=factor(), 
                            averages=double(), 
                            stDev=double(), 
                            tStats=double(), 
                            dfs=double(),
                            storedPs=double() )
  for ( i in 1:length( rois ) ) { # i <- 1
  
    dataSelectedLoop01 <- subset( dataSelected, roi==rois[i] & coefficient=='incongruent-congruent_GLT#0_Coef' )
    tLoop <- t.test( dataSelectedLoop01$medianValue )
    testOutputTemp <- data.frame( roi=dataSelectedLoop01$roi[1],
    				  roiIdx=dataSelectedLoop01$roiIdx[1],	
                                  averages=round( median( dataSelectedLoop01$medianValue ), 4 ), 
                                  stDev=round( sd( dataSelectedLoop01$medianValue ), 4 ), 
                                  tStats=tLoop$statistic, 
                                  dfs=tLoop$parameter,
                                  storedPs=tLoop$p.value )
    testOutput <- rbind( testOutput, testOutputTemp )
  }
  testOutput$storePsBonf <- p.adjust( testOutput$storedPs, method=c('bonferroni') )
  
  significantRois <- testOutput[ testOutput$storePsBonf < fdrThreshold, ]
  
  setwd( surfacesFolder )
  
  atlasFile <- read.AFNI( atlasSelected )
  atlasVolume <- atlasFile$brk[,,,2] # atlas information is stored in the second brik
  emptyVolume <- array(0,dim(atlasVolume))
  #roiNames <-  c('V1v',	'V1d', 'V2v', 'V2d', 'V3v', 'V3d', 'hV4', 'VO1', 'VO2', 'PHC1', 'PHC2',	    
  #    'MST', 'hMT', 'LO2', 'LO1', 'V3b', 'V3a', 'IPS0', 'IPS1', 'IPS2', 'IPS3', 'IPS4', 'IPS5', 'SPL1','FEF')	
  
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
  
  write.AFNI( filename = outputFile, 
              brk=emptyVolume,
              origin=atlasFile$origin,
              defhead=atlasFile$NI_head,
              orient=atlasFile$orient )
  
  return( significantRois )
  
}

#### results: output_ORIG_Denoised ####
setwd( mainFolder )
setwd( outputFolder )
sideSelected <- 'lh'
dataSelected <- subset( dataSelected_agg, side==sideSelected )
atlasSelected <- lhAtlas
outputFileLh <- 'output_Wang_lh_ORIG_Denoised.nii.gz'
fdrThreshold <- 0.05
sigRoiLh_orig_denoised <- generateFdrVolume( sideSelected, dataSelected, atlasSelected, outputFileLh, fdrThreshold, surfacesFolder )
print( sigRoiLh_orig_denoised )

setwd( mainFolder )
setwd( outputFolder )
sideSelected <- 'rh'
dataSelected <- subset( dataSelected_agg, side==sideSelected )
atlasSelected <- rhAtlas
outputFileRh <- 'output_Wang_rh_ORIG_Denoised.nii.gz'
fdrThreshold <- 0.05
sigRoiRh_orig_denoised <- generateFdrVolume( sideSelected, dataSelected, atlasSelected, outputFileRh, fdrThreshold, surfacesFolder )
print( sigRoiRh_orig_denoised )

setwd( surfacesFolder )
# clean up and combine output files
if ( file.exists( sprintf('output_Wang_ORIG_Denoised_combined.nii.gz',surfacesFolder) ) ) {
  instr <- sprintf( 'rm %s', sprintf('output_Wang_ORIG_Denoised_combined.nii.gz',surfacesFolder) )
  print( instr )
  system( instr )
}

setwd( surfacesFolder )
instr <- sprintf( '3dcalc -a %s -b %s -expr \u0027a+b\u0027 -prefix output_Wang_ORIG_Denoised_combined.nii.gz', outputFileLh, outputFileRh )
print( instr )
system( instr )


