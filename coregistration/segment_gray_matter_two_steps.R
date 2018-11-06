args <- commandArgs(T)
print( args )


# add a dilate_input step for gray matter only

#rm(list=ls())
#setwd('/data1/projects/myelin/myelinData/hemiBackup/hemianopticdata/V5980.hemianoptic')
#args <- c('pdCorrectRegularT1_noBlur_stripped.nii','white_matter_mask_0.09_0.01_10_7-7-3-3-5-5.nii.gz','/packages/afni/16.1.13','/data1/projects/myelin/analysisAfni/surfaces','/packages/afni/16.1.13','0.7-0.7-0.65-0.65')

#setwd('/data1/projects/myelin/myelinData/V2552')
#args <- c('pdCorrectRegularT1_noBlur_stripped.nii','white_matter_mask.nii.gz','/packages/afni/16.1.13','/data1/projects/myelin/analysisAfni/surfaces','/packages/afni/16.1.13','0.8','16')

#setwd('/data1/projects/myelin/myelinData/hemiBackup/hemianopticdata/V6435.hemianoptic')
#args <- c('pdCorrectRegularT1_noBlur_stripped_V6435_HE_straight.nii','rotated_white_matter_mask_0.01_0.09_10_5-5-2-2-0-1_0_6_1.nii.gz','/packages/afni/16.1.13','/data1/projects/myelin/analysisAfni/surfaces','/packages/afni/16.1.13','0.6-0.6-0.6-0.6')

source( sprintf('%s/AFNIio.R', args[3] ) )
source( sprintf('%s/coordinateFromLinearIndex.r', args[4] ) )
source( sprintf('%s/linearIndexFromCoordinate.r', args[4] ) )
library(RANN)
mainDir <- getwd()
atlasDir <- args[5]
imgInner <- args[2]
imgAmplitude <- args[1]

grayThrArray  <- as.numeric( strsplit( args[6], '[-]' )[[1]] ) #put here an array of 4 values
if (length(grayThrArray)!=4) { # missing argument
  msg <- sprintf( 'gray matter threshold should be 4 numbers' )
  warning( msg )
  stopifnot(flagDir)  
}

#nRoiClusters <- as.numeric( args[7] )

#print( minDiffThr )

instr <- sprintf('3dmask_tool -input %s -prefix dilate_wm.nii.gz -dilate_input 14', args[2] )
system( instr )

innerFile <- read.AFNI( imgInner )
innerVolume <- innerFile$brk[,,,1]
inside <- which( innerVolume==1 )
coordsInside <- coordinateFromLinearIndex( inside, dim(innerVolume) )

dilateFile <- read.AFNI( sprintf( 'dilate_wm.nii.gz' ) )
dilateVolume <- dilateFile$brk[,,,1]
dilateIndex <- which( dilateVolume == 1 & innerVolume == 0 )
coordsOutside <- coordinateFromLinearIndex( dilateIndex, dim(dilateVolume) )

iPart <- t( coordsInside )
oPart <- t( coordsOutside )
nnOut <- nn2( iPart, oPart, k=1 )
nnOutIn <- nn2( oPart, iPart, k=1 )

emptyVol <- array( 0, dim(innerVolume) )
emptyVol[inside] <- round( nnOutIn$nn.dists*-1, 2 )
emptyVol[dilateIndex] <- round( nnOut$nn.dists, 2 )
amplitudeData <- read.AFNI( imgAmplitude )
amplitudeData$brk[,,,1] <- emptyVol
volFileName <- sprintf( 'distanceMap.nii.gz' )
write.AFNI(volFileName,
           brk=emptyVol,
           label=NULL,
           view='+orig',
           orient=amplitudeData$orient,
           origin=amplitudeData$origin,
           defhead=amplitudeData$NI_head )

instr <- sprintf('conformToTlrc.sh %s %s/TT_N27+tlrc 0', args[1], atlasDir )
system( instr )
instr <- sprintf('cat_matvec MPRAGE_at.Xat.1D -I -ONELINE > invMat.1D' ) #invert the matrix
system( instr )

instr <- sprintf( '3dAutobox -prefix ttt_autoBox.nii.gz -noclust -input %s', imgAmplitude )
system( instr )

instr <- sprintf('3dAllineate -1Dmatrix_apply invMat.1D -source %s/TT_desai_dkpmaps+tlrc[0..15] -master ttt_autoBox.nii.gz -prefix atlas_map_01.nii.gz -final NN', atlasDir )
system( instr ) # get atlas part 01 in anatomical space
instr <- sprintf('3dAllineate -1Dmatrix_apply invMat.1D -source %s/TT_desai_dkpmaps+tlrc[16..32] -master ttt_autoBox.nii.gz -prefix atlas_map_02.nii.gz -final NN', atlasDir )
system( instr ) # get atlas part 02 in anatomical space
instr <- sprintf('3dAllineate -1Dmatrix_apply invMat.1D -source %s/TT_desai_dkpmaps+tlrc[33..48] -master ttt_autoBox.nii.gz -prefix atlas_map_03.nii.gz -final NN', atlasDir )
system( instr ) # get atlas part 03 in anatomical space
instr <- sprintf('3dAllineate -1Dmatrix_apply invMat.1D -source %s/TT_desai_dkpmaps+tlrc[49..69] -master ttt_autoBox.nii.gz -prefix atlas_map_04.nii.gz -final NN', atlasDir )
system( instr ) # get atlas part 04 in anatomical space

instr <- sprintf('3dmerge -1blur_fwhm 4.5 -doall -prefix atlas_map_01_blur.nii.gz atlas_map_01.nii.gz')
system( instr )
instr <- sprintf('3dmerge -1blur_fwhm 4.5 -doall -prefix atlas_map_02_blur.nii.gz atlas_map_02.nii.gz')
system( instr )
instr <- sprintf('3dmerge -1blur_fwhm 4.5 -doall -prefix atlas_map_03_blur.nii.gz atlas_map_03.nii.gz')
system( instr )
instr <- sprintf('3dmerge -1blur_fwhm 4.5 -doall -prefix atlas_map_04_blur.nii.gz atlas_map_04.nii.gz')
system( instr )

instr <- '3dcalc -a atlas_map_01_blur.nii.gz -expr \u0027step(a)\u0027 -prefix atlas_map_01_mask.nii.gz'
system( instr )
instr <- '3dcalc -a atlas_map_02_blur.nii.gz -expr \u0027step(a)\u0027 -prefix atlas_map_02_mask.nii.gz'
system( instr )
instr <- '3dcalc -a atlas_map_03_blur.nii.gz -expr \u0027step(a)\u0027 -prefix atlas_map_03_mask.nii.gz'
system( instr )
instr <- '3dcalc -a atlas_map_04_blur.nii.gz -expr \u0027step(a)\u0027 -prefix atlas_map_04_mask.nii.gz'
system( instr )

instr <- ('3dresample -master ttt_autoBox.nii.gz -prefix distance_map_res.nii.gz -rmode NN -input distanceMap.nii.gz')
system( instr )

#atlasMapName01 <- sprintf('atlas_map_01_mask.nii.gz') # load atlas ROIs (big files)
#atlasFile01 <- read.AFNI( atlasMapName01 )
#atlasVolume014D<- atlasFile01$brk

#atlasMapName02 <- sprintf('atlas_map_02_mask.nii.gz')
#atlasFile02 <- read.AFNI( atlasMapName02 )
#atlasVolume024D<- atlasFile02$brk 

#atlasMapName03 <- sprintf('atlas_map_03_mask.nii.gz')
#atlasFile03 <- read.AFNI( atlasMapName03 )
#atlasVolume034D<- atlasFile03$brk 

#atlasMapName04 <- sprintf('atlas_map_04_mask.nii.gz')
#atlasFile04 <- read.AFNI( atlasMapName04 )
#atlasVolume044D<- atlasFile04$brk 

innerFile <- read.AFNI( sprintf('distance_map_res.nii.gz') ) #load computed distance map
innerVolume <- innerFile$brk[,,,1]
innerVolume <- innerVolume
innerVolumeCopy <- innerVolume
innerVolumeCopy01 <- innerVolume

amplitudeData <- read.AFNI( 'ttt_autoBox.nii.gz' )
anatVolume <- amplitudeData$brk[,,,1] # load anatomical data

labIdxArray <- seq( 1 : 70 )

# areas lobe: 1=frontal, 2=temporal, 3=parietal, 4=occipital
lobes <- c(2,1,1,1,4,2,2,3,2,2,4,1,4,1,
           2,2,3,1,1,1,4,3,3,3,3,1,1,1,
           3,2,3,1,2,2,1,2,1,1,1,4,2,2,
           3,2,3,4,1,4,1,2,2,3,1,1,1,4,
           3,3,3,3,1,1,1,3,2,3,1,2,2,1)

emptyVolGray <- array( 0, dim(amplitudeData$brk[,,,1]) )
emptyVolLin <- array( 0, dim(amplitudeData$brk[,,,1]) )
emptyVolQuad <- array( 0, dim(amplitudeData$brk[,,,1]) )
#storeEstimatesCell <- list(0) #fix this

for ( ia in 1:length(labIdxArray) ) { # loop through atlas indexes, add control for numer of idx x ROI


  if (ia!=4 & ia!=39) { #do not analyze the corpus callosum
    
    print( sprintf('label: %d', ia) )

    if (ia == 1) {
      atlasMapName01 <- sprintf('atlas_map_01_mask.nii.gz') # load atlas ROIs (big files)
      atlasFile01 <- read.AFNI( atlasMapName01 )
      atlasVolume014D<- atlasFile01$brk
    }
    if (ia == 17) { 
      atlasMapName01 <- sprintf('atlas_map_02_mask.nii.gz') # load atlas ROIs (big files)
      atlasFile01 <- read.AFNI( atlasMapName01 )
      atlasVolume014D<- atlasFile01$brk
    }
    if (ia == 34) { 
      atlasMapName01 <- sprintf('atlas_map_03_mask.nii.gz') # load atlas ROIs (big files)
      atlasFile01 <- read.AFNI( atlasMapName01 )
      atlasVolume014D<- atlasFile01$brk
    }
    if (ia == 50) { 
      atlasMapName01 <- sprintf('atlas_map_04_mask.nii.gz') # load atlas ROIs (big files)
      atlasFile01 <- read.AFNI( atlasMapName01 )
      atlasVolume014D<- atlasFile01$brk
    }
    if (ia <= 16) { atlasVolume01 <- atlasVolume014D[,,,ia] }
    if (ia >= 17 & ia < 34) { atlasVolume01 <- atlasVolume014D[,,,ia-16] }
    if (ia >= 34 & ia < 50) { atlasVolume01 <- atlasVolume014D[,,,ia-33] }
    if (ia >= 50 ) { atlasVolume01 <- atlasVolume014D[,,,ia-49] }
    
    
    idx <- which( atlasVolume01 > 0.00001 & innerVolume !=0 & innerVolume > -1 & innerVolume < 8 & anatVolume!=0 )
    
    if (length(idx) > 30000) {
    
      nBinsClustering <- floor( length( idx ) / 28000 )
      nClusters <- nBinsClustering
      coordsIdx <- t( coordinateFromLinearIndex( idx, dim(anatVolume) ) )
      kIdx <- kmeans( coordsIdx , nClusters, iter.max = 450, algorithm="MacQueen" )
    
      storeGmEstimates <- rep(0,nClusters)
    
      for ( ik in 1:nClusters ) {
      
        idxK <- idx[ kIdx$cluster==ik ]
        anatIdx <- anatVolume[ idxK ]
        innerIdx <- innerVolume[ idxK ] + rnorm( length(idxK), 0, 0.05)
      
        nBinsQuantile <- floor( length( idxK ) / 750 )
        pQuantiles <- seq(0,1,length.out = nBinsQuantile)
        distBins <- .bincode( innerIdx, breaks=quantile( innerIdx, pQuantiles ), include.lowest = TRUE ) 
      

        # make it flexible with the number of voxels in the roi
        # bring out the computed quantile as a volume
        # post processing, compute distance again from gm mask
        # for each small subROI compute gm distribution, anf for the voxels
        # at the boundary [-1,+1] voxel, decide whether or not to add or eliminate them
        # cruise, do not give the hard segmentation but the probability (levelset to prob)
        # map quadratic and linear trend to the cortex

        distAvg <- tapply( innerIdx, list( distBins ), mean  ) 
        anatMean <- tapply( anatIdx, list( distBins ), mean  )
        anatSd01 <- tapply( anatIdx, list( distBins ), quantile, probs=c(0.025)  )
        anatSd02 <- tapply( anatIdx, list( distBins ), quantile, probs=c(0.975)  ) 
        anatSd <- anatSd02 - anatSd01
        #anatSd <- tapply( anatIdx, list( distBins ), sd  ) 
      
        startingPoint <- ceiling( length( unique( distBins ) ) / 2 )
        anatSdSmall <- anatSd[startingPoint:length(distAvg)]
        distAvgSmall <- distAvg[startingPoint:length(distAvg)]
        linPred <- seq(-1,1,length.out = length(anatSdSmall) )
        quadPred <- seq(-1,1,length.out = length(anatSdSmall) )^2 
        #qMode <- summary( lm( anatSdSmall ~ linPred + quadPred ) )
        #coeffLin <- coefficients( qMode )[2,3]
        #coeffQuad <- coefficients( qMode )[3,3]
        #pQuad <- coefficients( qMode )[3,4]
      
        ### run a second step using the quadratic index as a prior????
        qModeLin <- summary( lm( anatSdSmall ~ linPred ) )
        qModeQuad <- summary( lm( anatSdSmall ~ quadPred ) )
        coeffLin <- coefficients( qModeLin )[2,3]
        coeffQuad <- coefficients( qModeQuad )[2,3]
        pQuad <- coefficients( qModeQuad )[2,4]
      
        currGrayThr <- grayThrArray[ lobes[ ia ] ]
      
        if ( coeffLin <= 0 ) {
          storeIdx01 <- which( anatSdSmall < quantile( anatSdSmall, 1-currGrayThr ) )
        }
        if ( coeffLin > 0 ) {
          storeIdx01 <- which( anatSdSmall > quantile( anatSdSmall, currGrayThr ) )
        }
      
        estimatedGm <- distAvgSmall[ storeIdx01[1] ]
      
      
        #  if ( pQuad < 0.005  ) {
        #    estimatedGm <- distAvgSmall[ storeIdx01[1] ] + distAvgSmall[ storeIdx01[1] ]*0.05
        #  }
        #  if ( pQuad >= 0.005  ) {
        #    estimatedGm <- distAvgSmall[ storeIdx01[1] ]
        #  }
        
        gmValue <- array( estimatedGm )
        gmValue <- gmValue[1]
      
      
      
        print( round( gmValue, 2 ) )
        
        storeGmEstimates[ ik ] <- gmValue
      
        outVolValues <- emptyVolGray[ idxK ]
        gmIdx_zero <- idxK[ innerIdx >= 0 & innerIdx < gmValue & outVolValues==0 ]
        gmIdx <- idxK[ innerIdx >= 0 & innerIdx < gmValue & outVolValues!=0  ]
      
        if ( length( gmIdx_zero )>0 ) {
          emptyVolGray[ gmIdx_zero ] <- emptyVolGray[ gmIdx_zero ] + 0.1
          emptyVolLin[ gmIdx_zero ] <- emptyVolLin[ gmIdx_zero ] + coeffLin
          emptyVolQuad[ gmIdx_zero ] <- emptyVolQuad[ gmIdx_zero ] + coeffQuad
        }
        if ( length( gmIdx )>0 ) {
          emptyVolGray[ gmIdx ] <- ( emptyVolGray[ gmIdx ] + 0.1 ) / 2
          emptyVolLin[ gmIdx ] <- ( emptyVolLin[ gmIdx ] + coeffLin ) / 2
          emptyVolQuad[ gmIdx ] <- ( emptyVolQuad[ gmIdx ] + coeffQuad ) / 2
        }
      
      }
    
      #storeEstimatesCell[ia] <- storeGmEstimates
    
    }
  }
}

emptyVolWhole <- array( 0, c( dim(amplitudeData$brk[,,,1]), 3 ) )
emptyVolWhole[,,,1] <- emptyVolGray
emptyVolWhole[,,,2] <- emptyVolLin
emptyVolWhole[,,,3] <- emptyVolQuad
volFileName <- sprintf( 'gray_matter_mask_temp.nii.gz' )
write.AFNI(volFileName,
           brk=emptyVolWhole,
           label=NULL,
           view='+orig',
           orient=amplitudeData$orient,
           origin=amplitudeData$origin,
           defhead=amplitudeData$NI_head )

instr <- sprintf( '3dresample -master %s -input gray_matter_mask_temp.nii.gz -prefix gray_matter_mask.nii.gz -rmode NN', imgAmplitude )
system( instr )

system( sprintf( '3dcalc -a gray_matter_mask.nii.gz[0] -expr \u0027step(a)\u0027 -prefix gray_matter_mask_out.nii.gz' ) )

system( 'rm MPRAGE*' ) # clean atlas coregistration data
system( 'rm invMat.1D' ) # clean atlas coregistration data
system('rm pre.MPRAGE+orig*')
system('rm atlas_*')
system('rm dilate_wm.nii.gz')
system('rm distanceMap.nii.gz')
system('rm distance_map_res.nii.gz')
system('rm gray_matter_mask_temp.nii.gz')
system('rm gray_matter_mask.nii.gz')
system('rm ttt_autoBox.nii.gz')

