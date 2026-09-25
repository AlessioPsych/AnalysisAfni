rm( list=ls() )
source('~/abin/AFNIio.R')
library(lme4)
library(lmerTest)

# to run ' Rscript filename '


mainFolder <- '/mnt/disk01/ds000102_FlankerTask'
outputFolder <- 'derivatives/resultsGlasser'
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


lhAtlas <- 'lh.Glasser_HCP.nii.gz'
rhAtlas <- 'rh.Glasser_HCP.nii.gz'

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
  atlasVolume <- atlasFile$brk # atlas information is stored in the first brik
  emptyVolume <- array(0,dim(atlasVolume))
      #roiNames <-  c('V1','MST','V6','V2','V3','V4','V8','4','3b','FEF','PEF','55b','V3A','RSC','POS2','V7','IPS1','FFC','V3B','LO1','LO2','PIT','MT','A1','PSL','SFL',
      #               'PCV','STV','7Pm','7m','POS1','23d','v23ab','d23ab','31pv','5m','5mv','23c','5L','24dd','24dv','7AL','SCEF','6ma','7Am','7Pl','7PC','LIPv','VIP','MIP',
      #               '1','2','3a','6d','6mp','6v','p24pr','33pr','a24pr','p32pr','a24','d32','8BM','p32','10r','47m','8Av','8Ad','9m','8BL','9p','10d','8C','44','45','47l',
      #               'a47r','6r','IFJa','IFJp','IFSp','IFSa','p9-46v','46','a9-46v','9-46d','9a','10v','a10p','10pp','11l','13l','OFC','47s','LIPd','6a','i6-8','s6-8','43',
      #               'OP4','OP1','OP2-3','52','RI','PFcm','PoI2','TA2','FOP4','MI','Pir','AVI','AAIC','FOP1','FOP3','FOP2','PFt','AIP','EC','PreS','H','ProS','PeEc','STGa',
      #               'PBelt','A5','PHA1','PHA3','STSda','STSdp','STSvp','TGd','TE1a','TE1p','TE2a','TF','TE2p','PHT','PH','TPOJ1','TPOJ2','TPOJ3','DVT','PGp','IP2','IP1','IP0',
      #               'PFop','PF','PFm','PGi','PGs','V6A','VMV1','VMV3','PHA2','V4t','FST','V3CD','LO3','VMV2','31pd','31a','VVC','25','s32','pOFC','PoI1','Ig','FOP5','p10p','p47r',
      #               'TGv','MBelt','LBelt','A4','STSva','TE1m','PI','a32pr','p24')
  
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
outputFileLh <- 'output_Glasser_lh_ORIG_Denoised.nii.gz'
fdrThreshold <- 0.05
sigRoiLh_orig_denoised <- generateFdrVolume( sideSelected, dataSelected, atlasSelected, outputFileLh, fdrThreshold, surfacesFolder )
print( sigRoiLh_orig_denoised )

setwd( mainFolder )
setwd( outputFolder )
sideSelected <- 'rh'
dataSelected <- subset( dataSelected_agg, side==sideSelected )
atlasSelected <- rhAtlas
outputFileRh <- 'output_Glasser_rh_ORIG_Denoised.nii.gz'
fdrThreshold <- 0.05
sigRoiRh_orig_denoised <- generateFdrVolume( sideSelected, dataSelected, atlasSelected, outputFileRh, fdrThreshold, surfacesFolder )
print( sigRoiRh_orig_denoised )

setwd( surfacesFolder )
# clean up and combine output files
if ( file.exists( sprintf('output_Glasser_ORIG_Denoised_combined.nii.gz',surfacesFolder) ) ) {
  instr <- sprintf( 'rm %s', sprintf('output_Glasser_ORIG_Denoised_combined.nii.gz',surfacesFolder) )
  print( instr )
  system( instr )
}

setwd( surfacesFolder )
instr <- sprintf( '3dcalc -a %s -b %s -expr \u0027a+b\u0027 -prefix output_Glasser_ORIG_Denoised_combined.nii.gz', outputFileLh, outputFileRh )
print( instr )
system( instr )


