args <- commandArgs(T)
print( args )


# add a dilate_input step for gray matter only

#rm(list=ls())
#setwd('/home/fracasso/data/HighRes/Anatomies/Jelle')
#args <- c('anatomy_N3.nii.gz','white_matter_mask.nii.gz','/packages/afni/16.1.13','/home/fracasso/analysisAfni/surfaces','/packages/afni/16.1.13','-0.25')
args <- c('pdCorrectRegularT1_noBlur_stripped.nii','white_matter_mask_0.01_0.09_10_8-8-8-8-9-9.nii.gz','/packages/afni/16.1.13','/data1/projects/myelin/analysisAfni/surfaces','/packages/afni/16.1.13','-0.25')
setwd('/data1/projects/myelin/myelinData/hemiBackup')
setwd('./hemianopticdata/V6002.hemianoptic/SegmentGrayMatter_Test_JD/')

source( sprintf('%s/AFNIio.R', args[3] ) )
source( sprintf('%s/coordinateFromLinearIndex.r', args[4] ) )
source( sprintf('%s/linearIndexFromCoordinate.r', args[4] ) )
library(RANN)
mainDir <- getwd()
atlasDir <- args[5]
imgInner <- args[2]
imgAmplitude <- args[1]

minDiffThr <- args[6]

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
instr <- sprintf('3dAllineate -1Dmatrix_apply invMat.1D -source %s/TT_desai_dkpmaps+tlrc -master %s -prefix atlas_map_dkp.nii.gz -final NN', atlasDir, imgAmplitude )
system( instr ) # get atlas 02 in anatomical space

instr <- sprintf('3dTcat -prefix atlas_map_dkp_01.nii.gz  atlas_map_dkp.nii.gz[0..35]' )
system( instr )
instr <- sprintf('3dTcat -prefix atlas_map_dkp_02.nii.gz  atlas_map_dkp.nii.gz[36..69]' )
system( instr )

system( 'rm MPRAGE*' ) # clean atlas coregistration data
system( 'rm invMat.1D' ) # clean atlas coregistration data

atlasMapName01 <- sprintf('atlas_map_dkp_01.nii.gz') # load atlas ROIs (big files)
atlasFile01 <- read.AFNI( atlasMapName01 )
atlasVolume014D<- atlasFile01$brk

atlasMapName02 <- sprintf('atlas_map_dkp_02.nii.gz')
atlasFile02 <- read.AFNI( atlasMapName02 )
atlasVolume024D<- atlasFile02$brk 

innerFile <- read.AFNI( sprintf('distanceMap.nii.gz') ) #load computed distance map
innerVolume <- innerFile$brk[,,,1]
innerVolumeCopy <- innerVolume
innerVolumeCopy01 <- innerVolume

amplitudeData <- read.AFNI( imgAmplitude )
anatVolume <- amplitudeData$brk[,,,1] # load anatomical data

labIdxArray <- seq( 1 : c( dim(atlasVolume014D)[4]+dim(atlasVolume024D)[4] ) )

emptyVol <- array( 0, dim(amplitudeData$brk[,,,1]) )

for ( ia in 1:length(labIdxArray) ) { # loop through atlas indexes
#for ( ia in 1:36 ) { # loop through atlas indexes
  
  if (ia!=4 & ia!=39) { #do not analyze the corpus callosum
    
    print( sprintf('label: %d', ia) )
    if (ia <= 36) { atlasVolume01 <- atlasVolume014D[,,,ia] }
    if (ia > 36) { atlasVolume01 <- atlasVolume024D[,,,ia-36] }
    
    idx <- which( atlasVolume01 > 0.00001 & innerVolume !=0 & innerVolume > -1 & innerVolume < 7 & anatVolume!=0 )
    
    if ( length(idx)>4000 ) {
      nClusters <- 25
      coordsIdx <- t( coordinateFromLinearIndex( idx, dim(anatVolume) ) )
      kIdx <- kmeans( coordsIdx , nClusters, iter.max = 450, algorithm="MacQueen" )
    }
    if ( length(idx)<=4000 ) {
      nClusters <- 12
      coordsIdx <- t( coordinateFromLinearIndex( idx, dim(anatVolume) ) )
      kIdx <- kmeans( coordsIdx , nClusters, iter.max = 450, algorithm="MacQueen" )
    }
    
    #if ( min( table(kIdx$cluster) ) > 500 ) { 
    
    for ( ik in 1:nClusters ) {
      
      idxK <- idx[ kIdx$cluster==ik ]
      anatIdx <- anatVolume[ idxK ]
      innerIdx <- innerVolume[ idxK ]

      
      corticalDepthCorrelation <- function( xInput, yInput, amplitudeData ) {
        
        coordsDist <- coordinateFromLinearIndex( idxK, dim( emptyVol )  )
        cutBreakPoints <- quantile( xInput, seq(0,1,0.05) )
        depthCut <- as.numeric( .bincode( xInput, breaks = cutBreakPoints, include.lowest = TRUE ) )
        depthLevels <- tapply( xInput, list(depthCut), median  )
        sortedLevels <- sort( unique(depthCut) )
        storeCor <- rep(0, length(unique(depthCut))-1 )
        for (corCounter in 1:(length(unique(depthCut))-1) ) {
          coordDistSelStart <- t( coordsDist[, depthCut==sortedLevels[corCounter] ] )
          coordDistSelEnd <- t( coordsDist[, depthCut==sortedLevels[corCounter+1] ] )
          nnOut <- nn2( coordDistSelStart, coordDistSelEnd, k=1 )
        
          startCoords <- coordDistSelStart[nnOut$nn.idx,]
          endCoords <- coordDistSelEnd
          startIdx <- linearIndexFromCoordinate( t( startCoords ), dim( emptyVol )  )
          endIdx <- linearIndexFromCoordinate( t( endCoords ), dim( emptyVol )  )
          storeCor[corCounter] <- cor( amplitudeData$brk[ startIdx ], amplitudeData$brk[ endIdx ], method = c("pearson") )
        }
        
        return( cbind( storeCor, depthLevels[1:length(storeCor)] ) )
        #plot( storeCor ~ depthLevels[1:length(storeCor)] )
        
      }  
      
      corDF <- corticalDepthCorrelation( innerIdx, anatIdx, amplitudeData )
      storeIdx <- as.numeric( which( diff( corDF[,1] ) < minDiffThr ) )
      if ( length( storeIdx ) > 0 ) { 
        storeIdx2 <- which.max( storeIdx  )
        gmValue <- corDF[ storeIdx[storeIdx2] , 2 ] 
      }
      if ( length( storeIdx ) == 0 ) { 
        gmValue <- corDF[ dim(corDF)[1] , 2 ] 
      }
      
      print( gmValue )
      
      gmIdx <- idxK[ innerIdx >= 0 & innerIdx < gmValue ]
      emptyVol[ gmIdx ] <- emptyVol[ gmIdx ] + 0.1

      #gmIdx <- idxK[ innerIdx >= 0 & innerIdx < 5 ]
      #emptyVol[ gmIdx ] <- emptyVol[ gmIdx ] + 0.1
      
            
#        emptyVol[startIdx] <- 1
#        emptyVol[endIdx] <- 2
#        volFileName <- sprintf( 'delme.nii.gz' )
#        write.AFNI(volFileName,
#                   brk=emptyVol,
#                   label=NULL,
#                   view='+orig',
#                   orient=amplitudeData$orient,
#                   origin=amplitudeData$origin,
#                   defhead=amplitudeData$NI_head )
        
        
        
#        storeCor <- rep(0,9)
#        for ( idDepth in 1:9 ) {
#          storeCor[ idDepth ] <- cor( yInput[ depthCut==idDepth ], yInput[ depthCut==idDepth+1 ] )
#        } # try to match the numbers
#        return( storeCor )
        
#      }
      
      
      
      
      #sparArray <- seq(0.05,0.7,0.05)
      #smoothFun <- function( xInput, yInput, smoothInput ) {
      #  sp <- smooth.spline( xInput, yInput, spar=smoothInput )
      #  sDiff <- diff( sp$y )
        #ysDiff <- sp$y[2:length(sp$y)]
        #qThr <- quantile(yInput,0.5)
      #  indexGM01 <- which( sDiff == max( sDiff[ sDiff<=0 ] ) ) - 1
      #  gmValue <- sp$x[ indexGM01 ]
      #  sp$gmValue <- gmValue
      #  return( sp )
      #}
      #smoothOut <- lapply( sparArray, smoothFun, xInput=innerIdx, yInput=anatIdx)
      
      #indexPlot <- 5
      #plot( smoothOut[[indexPlot]]$y~smoothOut[[indexPlot]]$x, type='l' )
      #smoothOut[[indexPlot]]$gmValue
      
      #gmEst <- simplify2array( lapply( seq(1,length(sparArray)), function(dataIn,iList) { return(dataIn[[iList]]$gmValue) }, dataIn=smoothOut ) )
      
      ## from here! clean up the gmEst, they have to be all numerical
      
      #gmValue <- median( gmEst )
      #gmIdx <- idxK[ innerIdx >= 0 & innerIdx < gmValue ]
      #emptyVol[ gmIdx ] <- emptyVol[ gmIdx ] + 0.1
      #idxZero <- emptyVol[ gmIdx ] == 0
      #emptyVol[ gmIdx[ which( idxZero==TRUE ) ] ] <- 1 #FROM HERE!
      #emptyVol[ gmIdx[ which( idxZero==FALSE ) ] ] <- ( emptyVol[ gmIdx[ which( idxZero==FALSE ) ] ] + 1 ) / 2
      #sDiff <- diff( smoothOut[[indexPlot]]$y )
      #indexGM01 <- which( sDiff == max( sDiff[sDiff<0] ) ) - 1
      #gmValue <- sp$x[ indexGM01 ]
      #gmIdx <- idxK[ innerIdx > 0 & innerIdx <= gmValue ]
      #emptyVol[ gmIdx ] <- 1 
    }
  }
}

emptyVolOut <- array( 0, dim( emptyVol ) )
emptyVolOut[ emptyVol>0  ] <- 1
volFileName <- sprintf( 'gray_matter_mask.nii.gz' )
write.AFNI(volFileName,
           brk=emptyVolOut,
           label=NULL,
           view='+orig',
           orient=amplitudeData$orient,
           origin=amplitudeData$origin,
           defhead=amplitudeData$NI_head )

innerFile <- read.AFNI( imgInner )
whiteVol <- innerFile$brk[,,,1] 
emptyVolOut <- array( 0, dim( emptyVol ) )
emptyVolOut[ whiteVol>0 ] <- 2
emptyVolOut[ whiteVol==0 & emptyVol>0  ] <- 1
volFileName <- sprintf( 'segmentation_mask.nii.gz' )
write.AFNI(volFileName,
           brk=emptyVolOut,
           label=NULL,
           view='+orig',
           orient=amplitudeData$orient,
           origin=amplitudeData$origin,
           defhead=amplitudeData$NI_head )

system('rm pre.MPRAGE+orig*')
system('rm atlas_map_dkp.nii.gz')
system('rm atlas_map_dkp_01.nii.gz')
system('rm atlas_map_dkp_02.nii.gz')
system('rm dilate_wm.nii.gz')
system('rm distanceMap.nii.gz')
