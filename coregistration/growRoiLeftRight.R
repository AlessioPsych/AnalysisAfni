args <- commandArgs(T)
print( args )

setwd('/home/fracasso/data/HighRes/segmentationsAdjusted/segmentationAkhil')
args <- c('outputVolume.nii.gz','/usr/lib/afni/bin','/home/fracasso/analysisAfni/surfaces')

source( sprintf('%s/AFNIio.R', args[2] ) )
source( sprintf('%s/coordinateFromLinearIndex.r', args[3] ) )
source( sprintf('%s/linearIndexFromCoordinate.r', args[3] ) )
library(RANN)

roiFile <- read.AFNI(args[1]) 
roiVolume <- roiFile$brk[,,,1]
background <- which( roiVolume==0 )  
left <- which( roiVolume==1 )  
right <- which( roiVolume==2 )  

roiVolume[ left ] <- 0
roiFileName <- sprintf( '%s+orig', 'singleHemiRight' )
write.AFNI(roiFileName,
           brk=roiVolume,
           label=NULL,
           view='+orig',
           orient=roiFile$orient,
           origin=roiFile$origin,
           delta=roiFile$delta,
           defhead=roiFile$NI_head )
system('3dmerge -1blur_fwhm 5.0 -doall -prefix singleHemiRight_blur+orig singleHemiRight+orig')
roiFileBlur <- read.AFNI('singleHemiRight_blur+orig') 
volumeBlur <- roiFileBlur$brk[,,,1]
indxBlurRight <- volumeBlur > 0.1

roiFile <- read.AFNI(args[1]) 
roiVolume <- roiFile$brk[,,,1]
background <- which( roiVolume==0 )  
left <- which( roiVolume==1 )  
right <- which( roiVolume==2 )  

roiVolume[ right ] <- 0
roiVolume[ roiVolume!=0 ] <- 2
roiFileName <- sprintf( '%s+orig', 'singleHemiLeft' )
write.AFNI(roiFileName,
           brk=roiVolume,
           label=NULL,
           view='+orig',
           orient=roiFile$orient,
           origin=roiFile$origin,
           delta=roiFile$delta,
           defhead=roiFile$NI_head )
system( sprintf('3dmerge -1blur_fwhm %s.0 -doall -prefix singleHemiLeft_blur+orig singleHemiLeft+orig', args[4]) )
roiFileBlur <- read.AFNI('singleHemiLeft_blur+orig') 
volumeBlur <- roiFileBlur$brk[,,,1]
indxBlurLeft <- volumeBlur > 0.1

roiFile <- read.AFNI(args[1]) 
roiVolume <- roiFile$brk[,,,1]
emptyVolume <- array(0,dim(roiVolume))

emptyVolume[ indxBlurRight & roiVolume!=1 ] <- 1
emptyVolume[ indxBlurLeft & roiVolume!=2 ] <- 2

roiFileName <- 'left_right_smooth.nii.gz'
write.AFNI(roiFileName,
           brk=emptyVolume,
           label=NULL,
           view='+orig',
           orient=roiFile$orient,
           origin=roiFile$origin,
           delta=roiFile$delta,
           defhead=roiFile$NI_head )

system('rm singleHemiRight+orig.BRIK singleHemiRight+orig.HEAD singleHemiLeft+orig.BRIK singleHemiLeft+orig.HEAD')
system('rm singleHemiRight_blur+orig.BRIK singleHemiRight_blur+orig.HEAD singleHemiLeft_blur+orig.BRIK singleHemiLeft_blur+orig.HEAD')

# backgroundCoords <- t( coordinateFromLinearIndex( background, dim(roiVolume) ) )
# leftCoords <- t( coordinateFromLinearIndex( left, dim(roiVolume) ) )
# rightCoords <- t( coordinateFromLinearIndex( right, dim(roiVolume) ) )
# 
# nnOutLeft <- nn2( backgroundCoords,  leftCoords,  k=2 )
# nnOutRight <- nn2( backgroundCoords,  rightCoords,  k=1 )
# 
# selDistLeft <- nnOutLeft$nn.dists>=0 & nnOutLeft$nn.dists<=25
# selBackIndex <- nnOutLeft$nn.idx[selDistLeft]
# backIndex <- linearIndexFromCoordinate( t( backgroundCoords[selBackIndex,] ), dim( roiVolume ) )
# length(backIndex)
# 
# emptyVol <- array(0,dim(roiVolume))
# emptyVol[ backIndex ] <- 1
# 
# image( emptyVol[,,100] )
# 
# roiValue <- roiVolume[selBackIndex]



