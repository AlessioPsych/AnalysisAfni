args <- commandArgs(T)
print( args )

mainDir <- getwd()
setwd( mainDir )

#args <- c('V2551_HC_fitOutput.nii.gz',
#          'V2551_HC_pdCorrectRegularT1_noBlur_strippedcortexMask_calcarine_car_box_seg_depth.nii.gz',
#          'V2551_HC_fitOutput_expanded01.nii.gz',
#          '/packages/afni/17.0.13',
#          '/data1/projects/myelin/analysisAfni/surfaces')

source( sprintf('%s/AFNIio.R', args[4] ) )
source( sprintf('%s/coordinateFromLinearIndex.r', args[5] ) )
source( sprintf('%s/linearIndexFromCoordinate.r', args[5] ) )
library(RANN)

valuesBrk <- read.AFNI( args[1] )
depthBrk <- read.AFNI( args[2] )

valuesData <- valuesBrk$brk[,,,1]
depthData <- depthBrk$brk[,,,1]
volOut <- array( 0, c( dim( depthBrk$brk )[1:3], dim( valuesBrk$brk )[4] ) )

stepDepth <- 0.20
stepsFor <- seq(0,1,stepDepth)
for (nTimes in 1:(length( stepsFor )-1) ) {

  ind1 <- which( valuesData!=0 )
  ind2 <- which( depthData>=stepsFor[nTimes] & depthData<=stepsFor[nTimes+1] )
  coords01 <- matrix( t( coordinateFromLinearIndex( ind1, dim( depthData ) ) ), ncol=3 )
  coords02 <- matrix( t( coordinateFromLinearIndex( ind2, dim( depthData ) ) ), ncol=3 )
  nnOut <- nn2( coords02, coords01,  k=10 )
  
  valuesData <- valuesBrk$brk[,,,1]
  volOutTemp <- array( 0, c( dim( depthBrk$brk )[1:3] ) )
  
  valuesData <- valuesBrk$brk[,,,1]
  
  volOutTemp[ ind1 ] <- valuesData[ ind1 ]
  volOutTemp[ ind2[ nnOut$nn.idx[,1] ] ] <- valuesData[ ind1 ]
  volOutTemp[ ind2[ nnOut$nn.idx[,2] ] ] <- valuesData[ ind1 ]
  volOutTemp[ ind2[ nnOut$nn.idx[,3] ] ] <- valuesData[ ind1 ]
  volOutTemp[ ind2[ nnOut$nn.idx[,4] ] ] <- valuesData[ ind1 ]
  volOutTemp[ ind2[ nnOut$nn.idx[,5] ] ] <- valuesData[ ind1 ]
  #volOutTemp[ ind2[ nnOut$nn.idx[,6] ] ] <- valuesData[ ind1 ]
  #volOutTemp[ ind2[ nnOut$nn.idx[,7] ] ] <- valuesData[ ind1 ]
  #volOutTemp[ ind2[ nnOut$nn.idx[,8] ] ] <- valuesData[ ind1 ]
  #volOutTemp[ ind2[ nnOut$nn.idx[,9] ] ] <- valuesData[ ind1 ]
  #volOutTemp[ ind2[ nnOut$nn.idx[,10] ] ] <- valuesData[ ind1 ]
  
  if (nTimes==1) {
	volOut[,,,1] <- volOutTemp  
  }
  if (nTimes>1) {
	volOut[,,,1] <- volOut[,,,1] + volOutTemp  
  }
  
}
volOut[,,,1] <- volOut[,,,1]>0
idxClean <- depthData < 0
volOut[ idxClean ] <- 0

write.AFNI(args[3],
           brk=volOut,
           label=NULL,
           view='+orig',
           orient=valuesBrk$rient,
           origin=valuesBrk$origin,
           delta=valuesBrk$delta,
           defhead=valuesBrk$NI_head )




