args <- commandArgs(T)
print( args )

#setwd('/media/alessiofracasso/storage2/SpinozaTest/Numerosity_data_jelle/BaKl20160707')
#BOUNDARYNAME='boundary03'
#ANAT='anatomy_N3.nii.gz'
#DEPTH='depth.nii.gz'
#MAP='ModelingRes_NoZero1to7Log_al.nii'
#INDEXAVERAGE = '0'
#INDEXWEIGHT = '2'
#OUTNAME = 'name.nii.gz'
#args <- c( BOUNDARYNAME, ANAT, DEPTH, MAP, INDEXAVERAGE, INDEXWEIGHT, OUTNAME, '/usr/lib/afni/bin','/home/alessiofracasso/Dropbox/analysisAfni/surfaces')

source( sprintf('%s/AFNIio.R', args[8] ) )
source( sprintf('%s/coordinateFromLinearIndex.r', args[9] ) )
source( sprintf('%s/linearIndexFromCoordinate.r', args[9] ) )
library(RANN)
library(abind)

#get voxels on selected border
commandLine <- sprintf( '3dSurf2Vol -spec surfaces_folder/spec.surfaces.smoothed -surf_A surfaces_folder/%s_sm.1D.coord -sv %s -grid_parent %s -map_func mask -prefix %s', args[1], args[2], args[2], 'surfVolBoundary', args[1] )
print( commandLine )
system( commandLine )

borderMask <- read.AFNI('surfVolBoundary+orig')
depthMask <- read.AFNI( args[3] )

maskData <- borderMask$brk[,,,1]
depthData <- depthMask$brk[,,,1]

# connection matrix (storeIndx). Every voxel that intersects the relevant boundary, select voxels above and below (from WM to CSF)
# usign the depth map
stepLimit <- 0.1
limit1 <- seq(0,1-stepLimit,by=stepLimit)
limit2 <- limit1 + stepLimit
for (lim in 1:length(limit1) ) {
  ind1 <- which( maskData==1 )
  if (lim==1) {
    storeIndx <- array( 0, c( length(limit1), length(ind1) ) )
  }
  ind2 <- which( depthData>limit1[lim] & depthData<=limit2[lim] )
  coordsWM <- matrix( t( coordinateFromLinearIndex( ind1, dim(maskData) ) ), ncol=3 )
  coordsGM <- matrix( t( coordinateFromLinearIndex( ind2, dim(maskData) ) ), ncol=3 )
  nnOut <- nn2( coordsGM, coordsWM,  k=1 )
  coords <- coordsGM[ nnOut$nn.idx, ]
  coordsVol <- linearIndexFromCoordinate( t( coords ), dim( maskData ) )
  storeIndx[lim,] <- coordsVol
}

# resample map into parent space (probably the anatomy)
instr <- sprintf( '3dresample -master %s -inset %s -prefix upsampledMap.nii.gz -rmode NN', args[2], args[4] )
system( instr )

upsampledMap <- read.AFNI('upsampledMap.nii.gz')
upsampledMapData <- upsampledMap$brk
rm( list=c('upsampledMap','depthMask'))

# use the connection matrix (storeIndx) and extract the index and the weight matrix of neighboring voxels along the normals 
indexMat <- array( 0, dim( storeIndx ) )
indexWeight <- array( 0, dim( storeIndx ) )
matVolumeIndex <- as.numeric( args[5] ) + 1
matVolume <- upsampledMapData[,,,matVolumeIndex]
weightVolumeIndex <- as.numeric( args[6] ) + 1
weightVolume <- upsampledMapData[,,,weightVolumeIndex]
for (k in 1:dim(storeIndx)[1]) {
  indexMat[k,] <- matVolume[ storeIndx[k,] ] #index
  indexWeight[k,] <- weightVolume[ storeIndx[k,] ] #weight matrix
}

# compute the weighted index and store in into a volume (based on the parent data)
meanIndex <- apply( indexMat, 2, mean )
meanWeight <- apply( indexWeight, 2, mean )
wAverage <- apply( indexMat*indexWeight, 2, sum   ) / apply( indexWeight, 2, sum   )
ind1 <- which( maskData==1 )
emptyVol1 <- array( 0, dim(maskData) )
emptyVol2 <- array( 0, dim(maskData) )
emptyVol3 <- array( 0, dim(maskData) )
emptyVol1[ ind1 ] <- wAverage
emptyVol2[ ind1 ] <- meanIndex
emptyVol3[ ind1 ] <- meanWeight

roiFileName <- sprintf( '%s.nii.gz', 'wAverage' ) #save the weighted index volume
write.AFNI(roiFileName,
           brk=emptyVol1,
           label=NULL,
           view='+orig',           
           orient=borderMask$orient,
           origin=borderMask$origin,
           delta=borderMask$delta,
           defhead=borderMask$NI_head )

roiFileName <- sprintf( '%s.nii.gz', 'meanIndex' ) #save the mean index volume
write.AFNI(roiFileName,
           brk=emptyVol2,
           label=NULL,
           view='+orig',           
           orient=borderMask$orient,
           origin=borderMask$origin,
           delta=borderMask$delta,
           defhead=borderMask$NI_head )

roiFileName <- sprintf( '%s.nii.gz', 'meanWeight' ) #save the mean weight volume
write.AFNI(roiFileName,
           brk=emptyVol3,
           label=NULL,
           view='+orig',           
           orient=borderMask$orient,
           origin=borderMask$origin,
           delta=borderMask$delta,
           defhead=borderMask$NI_head )

emptyVol4 <- array( 0, dim(maskData) ) #save a volume to check how the connection matrix looks like
emptyVol4[ array( storeIndx[,seq(1,dim(storeIndx)[2],50)] ) ] <- 1
roiFileName <- sprintf( '%s.nii.gz', 'checkConnections' )
write.AFNI(roiFileName,
           brk=emptyVol4,
           label=NULL,
           view='+orig',           
           orient=borderMask$orient,
           origin=borderMask$origin,
           delta=borderMask$delta,
           defhead=borderMask$NI_head )

#concatenate the generated volumes with the upsampled map volume
instr <- sprintf( '3dTcat upsampledMap.nii.gz wAverage.nii.gz meanIndex.nii.gz meanWeight.nii.gz -prefix %s', args[7] )
system( instr )

#clean up a bit
system('rm surfVolBoundary+orig.BRIK')
system('rm surfVolBoundary+orig.HEAD')
system('rm wAverage.nii.gz')
system('rm meanIndex.nii.gz')
system('rm meanWeight.nii.gz')
system('rm upsampledMap.nii.gz')

