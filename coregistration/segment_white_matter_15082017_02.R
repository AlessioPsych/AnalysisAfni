args <- commandArgs(T)
print( args )

# free the iteration parameter

#rm(list=ls())
#setwd('/data1/projects/myelin/myelinData/hemiBackup/hemianopticdata/V5980.hemianoptic')
#args <- c('pdCorrectRegularT1_noBlur_stripped.nii','/packages/afni/17.0.13','/data1/projects/myelin/analysisAfni/surfaces','/packages/afni/17.0.13','0.001','0.07','5','6-6-6-6-6-6','/data1/projects/myelin/analysisAfni/coregistration','0','1','1')

# afni install dir / surfaces dir / atlas dir

source( sprintf('%s/AFNIio.R', args[2] ) )
source( sprintf('%s/coordinateFromLinearIndex.r', args[3] ) )
source( sprintf('%s/linearIndexFromCoordinate.r', args[3] ) )
source( sprintf('%s/computeDistanceMap.R',args[9]))
library(RANN)
library(ANTsR) # here is the ATNsR loading
mainDir <- getwd()
atlasDir <- args[4]
distanceThreshold <- as.numeric( strsplit( args[8], '[-]' )[[1]] ) #put here an array of 6 values

if (length(distanceThreshold)!=6) { # missing argument
  msg <- sprintf( 'distance array should be of 6 numbers' )
  warning( msg )
  stopifnot(flagDir)  
}
localProp <- as.numeric( args[5] )
globalProp <- as.numeric( args[6] )
biasVoxelDistance <- as.numeric( args[10] )
parametersInput <- as.numeric( args[11] )
paramsCleanUp <- as.numeric( args[12] )

counterManualStart <- as.numeric( strsplit( args[11], '[-]' )[[1]] ) 
classesManualStart <- as.numeric( strsplit( args[12], '[-]' )[[1]] ) 
counterManual <- rep(0,6)
counterManual[ counterManualStart ] <- classesManualStart

system( sprintf('cp %s/TT_N27_EZ_LR+tlrc.BRIK.gz %s', atlasDir, mainDir ) )
system( sprintf('cp %s/TT_N27_EZ_LR+tlrc.HEAD %s', atlasDir, mainDir ) )
system( sprintf('cp %s/TT_desai_dkpmaps+tlrc.HEAD %s', atlasDir, mainDir ) )
system( sprintf('cp %s/TT_desai_dkpmaps+tlrc.BRIK.gz %s', atlasDir, mainDir ) )
system( sprintf('cp %s/TT_desai_dkpmaps+tlrc.HEAD %s', atlasDir, mainDir ) )
system( sprintf('cp %s/TT_desai_dk_mpm+tlrc.BRIK.gz %s', atlasDir, mainDir ) )
system( sprintf('cp %s/TT_desai_dk_mpm+tlrc.HEAD %s', atlasDir, mainDir ) )
system( '3dTstat -sum -prefix atlasSum+tlrc TT_desai_dkpmaps+tlrc' )
system( '3dcalc -a atlasSum+tlrc -b TT_desai_dk_mpm+tlrc -expr \u0027 step( step(a) + WITHIN(b,89,89) + WITHIN(b,71,71) ) \u0027 -prefix atlasSumMask+tlrc')
system( '3dmask_tool -input atlasSumMask+tlrc -prefix atlasSumMaskFilled+tlrc -fill_holes' )

instr <- sprintf('conformToTlrc.sh %s %s/TT_N27+tlrc 0', args[1], atlasDir )
system( instr )
instr <- sprintf('cat_matvec MPRAGE_at.Xat.1D -I -ONELINE > invMat.1D' ) #invert the matrix
system( instr )
instr <- sprintf('3dAllineate -1Dmatrix_apply invMat.1D -source atlasSumMaskFilled+tlrc -master %s -prefix %s -final NN', args[1], 'mask.nii.gz' )
system( instr ) # mask brain only

system( 'rm MPRAGE*' ) # clean atlas coregistration data
system( 'rm invMat.1D' ) # clean atlas coregistration data
system( 'rm TT_N27**' )
system( 'rm TT_desai**' )

#instr <- sprintf('conformToTlrc.sh %s %s/MNI_avg152T1+tlrc 1', args[1], atlasDir )
#system( instr )
#instr <- sprintf('3dcalc -a MPRAGE_at.nii -b %s/MNI_avg152T1+tlrc -expr \u0027step(a)*within(x,0,1000)*1+step(a)*within(x,-1000,-0.001)*2\u0027 -prefix lr.nii.gz', atlasDir)
#system( instr )
#instr <- sprintf('convertRoi_from_desi_atlas.sh MPRAGE_at.Xat.1D %s lr.nii.gz 1 ', args[1])
#system(instr)

system( 'rm MPRAGE*' ) # clean atlas coregistration data
system( 'rm pre*' )

instr <- sprintf('3dmask_tool -input mask.nii.gz -prefix mask_dil.nii.gz -dilate_input 3')
system( instr )

### mask input to esclude cerebellum
instr <- sprintf('3dcalc -a %s -b mask_dil.nii.gz -expr \u0027step(b)*a\u0027 -prefix inputAnat_mask.nii.gz', args[1] )
system( instr )
### autobox volume
instr <- '3dAutobox -noclust -prefix inputAnat_mask_box.nii.gz -input inputAnat_mask.nii.gz'
system( instr )

# find the centers, the space and set the limits accordingly
system( sprintf( '3dinfo -o3 inputAnat_mask_box.nii.gz > origin.1D' ) )
system( sprintf( '3dinfo -n4 inputAnat_mask_box.nii.gz > dimensions.1D' ) )
system( sprintf( '3dinfo -ad3 inputAnat_mask_box.nii.gz > voxSize.1D' ) )
system( sprintf( '3dinfo -orient inputAnat_mask_box.nii.gz > orient.1D' ) )
origin <- as.numeric( scan( file= 'origin.1D', skip=0, nlines=1, what='character' ) )
dimensions4D <- as.numeric( scan( file= 'dimensions.1D', skip=0, nlines=1, what='character' ) )
dimensions3D <- dimensions4D[1:3]
voxsize <- as.numeric( scan( file= 'voxSize.1D', skip=0, nlines=1, what='character' ) )
orientation <- scan( file= 'orient.1D', skip=0, nlines=1, what='character' )
orientationSplit <- strsplit( orientation, split='' )[[1]]
volLimit <- origin + voxsize*dimensions3D 
index01 <- which( orientationSplit == 'A')
index02 <- which( orientationSplit == 'P')
index03 <- which( orientationSplit == 'I')
index04 <- which( orientationSplit == 'S')

if ( length( index01 )!=0 ) {
  antLimit <- origin[ index01 ]
  posLimit <- volLimit[ index01 ]
}
if ( length( index02 )!=0 ) {
  antLimit <- volLimit[ index02 ]
  posLimit <- origin[ index02 ]
}
if ( length( index03 )!=0 ) {
  infLimit <- origin[ index03 ]
  supLimit <- volLimit[ index03 ]
}
if ( length( index04 )!=0 ) {
  infLimit <- volLimit[ index04 ]
  supLimit <- origin[ index04 ]
}

system('rm orient.1D origin.1D dimensions.1D voxSize.1D')

limits <- round( seq( posLimit, antLimit, length.out = 7  ), 0 )
limitsIS <- round( seq( infLimit, supLimit, length.out = 3  ), 0 )

if (parametersInput==1) {
  parameterSegmentation <- read.table( 'paramFileAnatomy.1D', as.is=TRUE )
  parameterSegmentationProp <- parameterSegmentation[1:3,]
  parameterSegmentationClass <- parameterSegmentation[4:6,]
}
if (parametersInput==0) {
  parameterSegmentationProp <- array(0,c(3,6))
  parameterSegmentationClass <- array(3,c(3,6))
}


for ( counter in 1:(length(limits)-1) ) {
  #for ( counter in 1:2 ) {
  setwd(mainDir)
  if (counter==1) {
    system( sprintf( '@clip_volume -input %s -anterior %1.0f -prefix clip_0%1.0f.nii.gz', 'inputAnat_mask_box.nii.gz', limits[2], counter ) )
  }
  if (counter>1) {
    system( sprintf( '@clip_volume -input %s -posterior %1.0f -anterior %1.0f -prefix clip_0%1.0f.nii.gz', 'inputAnat_mask_box.nii.gz', limits[counter], limits[counter+1], counter ) )
  }
  instr <- sprintf('3dcalc -a clip_0%1.0f.nii.gz -expr \u0027step(a)\u0027 -prefix maskVolume.nii.gz', counter)
  system( instr )
  
  
  clipVolumeLoad <- read.AFNI( sprintf( 'clip_0%1.0f.nii.gz', counter ) ) # this is the quantile clipping
  clipVolumeLoadData <- clipVolumeLoad$brk[,,,1]
  thrClip <- quantile( array( clipVolumeLoadData ), probs=c(0.0001,0.9999) )
  clipVolumeLoadDataLogical <- clipVolumeLoadData>thrClip[1] &  clipVolumeLoadData<thrClip[2]
  clipVolumeLoadDataThr <- array( 0, dim(clipVolumeLoadData) )
  clipVolumeLoadDataThr[ clipVolumeLoadDataLogical ] <- clipVolumeLoadData[ clipVolumeLoadDataLogical ]
  clipVolumeLoad$brk[,,,1] <- clipVolumeLoadDataThr
  write.AFNI( filename=sprintf( 'clip_0%1.0f_thr.nii.gz', counter ),
              brk = clipVolumeLoad$brk,
              origin = clipVolumeLoad$origin, 
              delta = clipVolumeLoad$delta, 
              defhead = clipVolumeLoad$NI_head )
  
  volumeAnts <- antsImageRead( sprintf( 'clip_0%1.0f_thr.nii.gz', counter ) ); #here I use ANTsR, modify
  volumeAntsN3 <- n3BiasFieldCorrection( volumeAnts, 5 )
  antsImageWrite( volumeAntsN3, 'volumeN3.nii.gz' ) 
  
  instr <- sprintf('3dSeg -anat volumeN3.nii.gz -mask maskVolume.nii.gz -classes \u0027CSF ; GM ; WM\u0027 -bias_classes \u0027GM ; WM\u0027 -bias_fwhm 0.0 -mixfrac UNI -main_N %s -blur_meth BFT', args[7] ) 
  print( instr )
  system( instr )
  system('rm maskVolume.nii.gz')
  system( sprintf('mv Segsy Segsy_0%1.0f', counter) )
  
  instr <- sprintf('rm clip_0%1.0f.nii.gz', counter)
  system( instr )
  instr <- sprintf('rm clip_0%1.0f_thr.nii.gz', counter)
  system( instr )
  
  system('rm volumeN3.nii.gz')
  
  setwd( sprintf('Segsy_0%1.0f', counter) )
  
  # anatVol <- read.AFNI( 'Anat+orig' )
  # anatVol3d <- anatVol$brk[,,,1]
  # classVol <- read.AFNI( 'Classes+orig' )
  # classData3d <- classVol$brk[,,,1]
  # uniqueClasses <- unique( array( classData3d ) )
  # classMedian <- function(anatIn, volumeIn, classId) {
  #   idxVol <- volumeIn == classId
  #   out <- c( median( anatIn[ idxVol ] ), length( anatIn[ idxVol ]  ), classId )
  #   return( out )
  # }
  # classesMedian <- sapply( uniqueClasses, classMedian, anatIn = anatVol3d, volumeIn = classData3d )
  # classesProps <- prop.table( classesMedian[2,] ) * 100
  # classesMedianFilt <- classesMedian[ , classesProps > 0.5 ]
  # idxMaxClass <- which.max( classesMedianFilt[1,] )
  # 
  # if (counterManual[counter]==0) { classWM <- classesMedianFilt[ 3, idxMaxClass ] }
  # if (counterManual[counter]!=0) { classWM <- counterManual[counter] }
  
  
  
  
  
  anatVol <- read.AFNI( 'Anat+orig' )
  anatVol3d <- anatVol$brk[,,,1]
  classVol <- read.AFNI( 'Classes+orig' )
  classData3d <- classVol$brk[,,,1]
  
  instr <- '3dAutobox -input Posterior+orig -prefix PosteriorAuto.nii.gz -noclust'
  system( instr )
  instr <- '3dresample -master PosteriorAuto.nii.gz -input Anat+orig -prefix anatAuto.nii.gz -rmode NN'
  system( instr )
  instr <- '3dresample -master PosteriorAuto.nii.gz -input Classes+orig -prefix classesAuto.nii.gz -rmode NN'
  system( instr )
  
  system( sprintf( '3dinfo -o3 anatAuto.nii.gz > origin.1D' ) )
  system( sprintf( '3dinfo -n4 anatAuto.nii.gz > dimensions.1D' ) )
  system( sprintf( '3dinfo -ad3 anatAuto.nii.gz > voxSize.1D' ) )
  system( sprintf( '3dinfo -orient anatAuto.nii.gz > orient.1D' ) )
  origin <- as.numeric( scan( file= 'origin.1D', skip=0, nlines=1, what='character' ) )
  dimensions4D <- as.numeric( scan( file= 'dimensions.1D', skip=0, nlines=1, what='character' ) )
  dimensions3D <- dimensions4D[1:3]
  voxsize <- as.numeric( scan( file= 'voxSize.1D', skip=0, nlines=1, what='character' ) )
  orientation <- scan( file= 'orient.1D', skip=0, nlines=1, what='character' )
  orientationSplit <- strsplit( orientation, split='' )[[1]]
  volLimit <- origin + voxsize*dimensions3D 
  index01 <- which( orientationSplit == 'A')
  index02 <- which( orientationSplit == 'P')
  index03 <- which( orientationSplit == 'I')
  index04 <- which( orientationSplit == 'S')
  
  if ( length( index01 )!=0 ) {
    antLimitPart <- origin[ index01 ]
    posLimitPart <- volLimit[ index01 ]
  }
  if ( length( index02 )!=0 ) {
    antLimitPart <- volLimit[ index02 ]
    posLimitPart <- origin[ index02 ]
  }
  if ( length( index03 )!=0 ) {
    infLimitPart <- origin[ index03 ]
    supLimitPart <- volLimit[ index03 ]
  }
  if ( length( index04 )!=0 ) {
    infLimitPart <- volLimit[ index04 ]
    supLimitPart <- origin[ index04 ]
  }
  
  wmStep01 <- read.AFNI( 'PosteriorAuto.nii.gz' )
  wmStep01Volume <- wmStep01$brk[,,,1]
  limitsSize <- round( seq(infLimitPart,supLimitPart,length.out = 4), 0 )
  emptyVol <- array(0, dim( wmStep01Volume ) )
  outVol <- array(0, dim( wmStep01Volume ) )
  for (kSplit in 1:(length( limitsSize )-1) ) {
    classWM <- parameterSegmentationClass[kSplit,counter]
    if (parameterSegmentationProp[kSplit,counter]!=0) {
      system( sprintf( '@clip_volume -input %s -inferior %1.0f -superior %1.0f -prefix clip_0%1.0f.nii.gz', 'PosteriorAuto.nii.gz', limitsSize[kSplit], limitsSize[kSplit+1], kSplit ) )
      tempFile <- read.AFNI( sprintf('clip_0%1.0f.nii.gz', kSplit) )
      emptyVol <- tempFile$brk[,,,classWM]
      outVol[ emptyVol > parameterSegmentationProp[kSplit,counter]  ] <- 1
    }
    if (parameterSegmentationProp[kSplit,counter]==0) {
      system( sprintf( '@clip_volume -input %s -inferior %1.0f -superior %1.0f -prefix clip_0%1.0f.nii.gz', 'classesAuto.nii.gz', limitsSize[kSplit], limitsSize[kSplit+1], kSplit ) )
      tempFile <- read.AFNI( sprintf('clip_0%1.0f.nii.gz', kSplit) )
      emptyVol <- tempFile$brk[,,,1]
      outVol[ emptyVol == classWM  ] <- 1
    }
  }
  
  
  
  
  
  
  #emptyVol <- outVol
  #classVol$brk[,,,1] <- emptyVol
  write.AFNI( filename='wmMask.nii.gz',
              brk = outVol,
              origin = wmStep01$origin, 
              delta = wmStep01$delta, 
              defhead = wmStep01$NI_head )
  # isolate WM:
  #instr <- '3dcalc -a Classes+orig -expr \u0027within(a,3,3)\u0027 -prefix wmMask.nii.gz'
  #system( instr )
  instr <- '3dclust 0 20 wmMask.nii.gz > out.1D'
  system( instr )
  clustTable <- read.table( 'out.1D', comment.char = "#" )
  system( 'rm out.1D' )
  maxClust <- clustTable[1,1]
  clustProp <- round( clustTable[,1] / maxClust, 5 )
  whichIdxArray <- which( clustProp > localProp )
  whichIdx <- whichIdxArray[ length( whichIdxArray ) ]
  clusterCutOff <- clustTable[whichIdx,1]
  clusteringInst <- sprintf('3dclust -prefix wmMask_clust.nii.gz 0 %1.0f wmMask.nii.gz', clusterCutOff - 1 )
  system( clusteringInst )
  
  ### do it for the global step, check the results....
  
  system( '3dmask_tool -input wmMask_clust.nii.gz -fill_holes -prefix wmMask_clust_fill.nii.gz' )
  
  system( '3dmask_tool -input wmMask_clust_fill.nii.gz -dilate_input -1 -prefix wmMask_clust_fill01.nii.gz' ) #erosion
  
  
  imgOriginal <- 'wmMask_clust_fill.nii.gz'
  innerOriginal <- read.AFNI( imgOriginal )
  innerVolumeOriginal <- innerOriginal$brk[,,,1]
  
  imgInner <- 'wmMask_clust_fill01.nii.gz'
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
  volFileName <- sprintf( 'distanceMap.nii.gz' )
  write.AFNI(volFileName,
             brk=emptyVol,
             label=NULL,
             view='+orig',
             orient=innerFile$orient,
             origin=innerFile$origin,
             defhead=innerFile$NI_head )
  
  indexDeleteFromOriginal <- which( emptyVol>distanceThreshold[ counter ] ) #if distance larger that threshold, then set to zero
  innerVolumeOriginal[ indexDeleteFromOriginal ] <- 0
  innerOriginal$brk[,,,1] <- innerVolumeOriginal
  volFileName <- sprintf( 'wmMask_clust_fill_distance.nii.gz' )
  write.AFNI(volFileName,
             brk=innerOriginal$brk,
             label=NULL,
             view='+orig',
             orient=innerOriginal$orient,
             origin=innerOriginal$origin,
             defhead=innerOriginal$NI_head )
  
  instr <- '3dclust 0 20 wmMask_clust_fill_distance.nii.gz > out.1D'
  system( instr )
  clustTable <- read.table( 'out.1D', comment.char = "#" )
  system( 'rm out.1D' )
  maxClust <- clustTable[1,1]
  clustProp <- round( clustTable[,1] / maxClust, 5 )
  whichIdxArray <- which( clustProp > localProp )
  whichIdx <- whichIdxArray[ length( whichIdxArray ) ]
  clusterCutOff <- clustTable[whichIdx,1]
  clusteringInst <- sprintf('3dclust -prefix wmMask_clust_fill_distance_clust.nii.gz 0 %1.0f wmMask_clust_fill_distance.nii.gz', clusterCutOff - 1 )
  system( clusteringInst )
  
  system( '3dmask_tool -input wmMask_clust_fill_distance_clust.nii.gz -fill_holes -prefix wmMask_clust_fill_distance_clust_fill.nii.gz' )
  
  system( '3dmask_tool -input wmMask_clust_fill_distance_clust_fill.nii.gz -dilate_input 1 -1 -prefix wmMask_clust_fill_distance_clust_fill_dil_er.nii.gz' ) #erosion
  
  system( '3dresample -master Anat+orig -input wmMask_clust_fill_distance_clust_fill_dil_er.nii.gz -prefix wmMask_clust_fill_distance_clust_fill_dil_er01.nii.gz' ) #erosion
  
  #source( sprintf('%s/computeDistanceMap.R',args[9]))
  
  #computeDistanceMap( 'wmMask_clust_fill_distance_clust_fill_dil_er.nii.gz' )
  
}

setwd(mainDir)
dirSegmentation <- dir( pattern=sprintf( 'Segsy*') )
letterArray <- letters[1:length( dirSegmentation )]
for (k in 1:length( dirSegmentation ) ) {
  if (k==1) { instrVar <- sprintf( '-%s %s/wmMask_clust_fill_distance_clust_fill_dil_er01.nii.gz', letterArray[k], dirSegmentation[k] ) }
  if (k>1) { instrVar <- sprintf( '%s -%s %s/wmMask_clust_fill_distance_clust_fill_dil_er01.nii.gz', instrVar, letterArray[k], dirSegmentation[k] ) }
  if (k==1) { instrExpr <- sprintf( '\u0027%s', letterArray[k] ) }
  if (k>1) { instrExpr <- sprintf( '%s + %s',instrExpr, letterArray[k] ) } 
  if (k==length( dirSegmentation )) { instrExpr <- sprintf( '%s\u0027', instrExpr, letterArray[k] ) } 
}
instrCalc <- sprintf('3dcalc %s -expr %s -prefix wm_segVol.nii.gz', instrVar, instrExpr)
system( instrCalc )

# isolate WM:
instr <- '3dclust 0 50 wm_segVol.nii.gz > out.1D'
system( instr )
clustTable <- read.table( 'out.1D', comment.char = "#" )
system( 'rm out.1D' )
maxClust <- clustTable[1,1]
clustProp <- round( clustTable[,1] / maxClust, 5 )
whichIdxArray <- which( clustProp > globalProp )
whichIdx <- whichIdxArray[ length( whichIdxArray ) ]
clusterCutOff <- clustTable[whichIdx,1]
clusteringInst <- sprintf('3dclust -prefix wmMask_clust.nii.gz 0 %1.0f wm_segVol.nii.gz', clusterCutOff - 1 )
system( clusteringInst )

system( '3dmask_tool -input wmMask_clust.nii.gz -fill_holes -prefix wmMask_clust_fill.nii.gz' )
system( '3dmask_tool -input wmMask_clust_fill.nii.gz -dilate_input 1 -1 -prefix wmMask_clust_fill01.nii.gz' )
system( '3dmask_tool -input wmMask_clust_fill01.nii.gz -fill_holes -prefix wmMask_clust_fill02.nii.gz' )

computeDistanceMap('wmMask_clust_fill02.nii.gz')

if (biasVoxelDistance == 0) {
  system('3dcopy wmMask_clust_fill02.nii.gz wmMask_clust_fill03.nii.gz')
}
if (biasVoxelDistance > 0) {
  biasDistance <- abs( as.numeric( biasVoxelDistance  ) )
  instr <- sprintf('3dcalc -a ttt_distanceMap_overall+orig -expr \u0027within(a-%s,-1000,0)*step(a+10000)\u0027 -prefix wmMask_clust_fill03.nii.gz', biasDistance)
  system( instr )
}
if (biasVoxelDistance < 0) {
  biasDistance <- abs( as.numeric( biasVoxelDistance  ) )
  instr <- sprintf('3dcalc -a ttt_distanceMap_overall+orig -expr \u0027within(a+%s,-1000,0)*step(a+10000)\u0027 -prefix wmMask_clust_fill03.nii.gz', biasDistance)
  system( instr )
}

system( sprintf( '3dresample -inset wmMask_clust_fill03.nii.gz -prefix wmMask_clust_fill04.nii.gz -master %s -rmode NN', args[1] ) )

if (paramsCleanUp==1) {
  system('rm atlas*')
  system('rm mask*')
  system('rm wmMask_clust.nii.gz')
  system('rm wm_segVol.nii.gz')
  system('rm -R Segsy*')
  system('rm inputAnat*')
  system('rm ttt_*')
  system('rm wmMask_clust_fill.nii.gz')
  system('rm wmMask_clust_fill01.nii.gz')
  system('rm wmMask_clust_fill02.nii.gz')
  system('rm wmMask_clust_fill03.nii.gz')
  #system( sprintf( '3dresample -master %s -inset wmMask_clust_fill03.nii.gz -prefix white_matter_mask.nii.gz',args[1] ) ) 
  system('mv wmMask_clust_fill04.nii.gz white_matter_mask.nii.gz')
}