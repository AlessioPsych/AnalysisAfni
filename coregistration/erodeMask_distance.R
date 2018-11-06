args <- commandArgs(T)
print( args )

#args <- c('braimask_corr.nii.gz','3','delMe.nii.gz')

source( sprintf('%s/AFNIio.R', Sys.getenv(x='AFNI_INSTALLDIR') ) )
source( sprintf('%s/coordinateFromLinearIndex.r', Sys.getenv(x='AFNI_TOOLBOXDIRSURFACES') ) )
source( sprintf('%s/linearIndexFromCoordinate.r', Sys.getenv(x='AFNI_TOOLBOXDIRSURFACES') ) )
source( sprintf('%s/computeDistanceMap.R', Sys.getenv(x='AFNI_TOOLBOXDIRCOREGISTRATION')  ))
library(RANN)

inputMask <- args[1]
distanceThr <- as.numeric( args[2] )
outputFileName <- args[3]

instr <- sprintf( '3dcopy %s ttt_mask_copy00.nii.gz', inputMask)
system( instr )
system( '3dmask_tool -input ttt_mask_copy00.nii.gz -dilate_input -1 -prefix ttt_mask_copy01.nii.gz' ) #erosion

imgOriginal <- 'ttt_mask_copy00.nii.gz'
innerOriginal <- read.AFNI( imgOriginal )
innerVolumeOriginal <- innerOriginal$brk[,,,1]

imgInner <- 'ttt_mask_copy01.nii.gz'
innerFile <- read.AFNI( imgInner )
innerVolume <- innerFile$brk[,,,1]
inside <- which( innerVolume==1 )

voxelDel <- innerVolume != innerVolumeOriginal
outside <- which( voxelDel==1 ) 

coordsInside <- coordinateFromLinearIndex( inside, dim(innerVolume) )
coordsOutside <- coordinateFromLinearIndex( outside, dim(innerVolume) )

iPart <- t( coordsInside ) # distances from eroded volume to eroded voxels
oPart <- t( coordsOutside ) 
nnOut <- nn2( iPart, oPart, k=1 )
nnOutIn <- nn2( oPart, iPart, k=1 )

emptyVol <- array( 0, dim(innerVolume) )
emptyVol[inside] <- round( nnOutIn$nn.dists*-1, 2 )
emptyVol[outside] <- round( nnOut$nn.dists, 2 )
innerFile$brk[,,,1] <- emptyVol
volFileName <- sprintf( 'ttt_distanceMap.nii.gz' )
write.AFNI(volFileName,
           brk=emptyVol,
           label=NULL,
           view='+orig',
           orient=innerFile$orient,
           origin=innerFile$origin,
           defhead=innerFile$NI_head )

indexDeleteFromOriginal <- which( emptyVol>distanceThr ) #if distance larger that threshold, then set to zero
innerVolumeOriginal[ indexDeleteFromOriginal ] <- 0
innerOriginal$brk[,,,1] <- innerVolumeOriginal
volFileName <- sprintf( 'ttt_mask_copy02.nii.gz' )
write.AFNI(volFileName,
           brk=innerOriginal$brk,
           label=NULL,
           view='+orig',
           orient=innerOriginal$orient,
           origin=innerOriginal$origin,
           defhead=innerOriginal$NI_head )

system( sprintf( 'mv ttt_mask_copy02.nii.gz %s', outputFileName ) )
system('rm ttt_mask_copy00.nii.gz')
system('rm ttt_mask_copy01.nii.gz')
system('rm ttt_distanceMap.nii.gz')
