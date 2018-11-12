#args <- commandArgs(T)
#print( args )


rm(list=ls())
setwd('/data1/projects/myelin/myelinData/hemiBackup/hemianopticdata/V6512.hemianoptic')
args <- c('pdCorrectRegularT1_noBlur_stripped.nii',
          '/packages/afni/17.0.13','/data1/projects/myelin/analysisAfni/surfaces','/packages/afni/17.0.13',
          '6.5-6.5-7-7-6.5-6.5','10','white_matter_mask_0.01_0.09_20_5-5-5-8-8-8.nii.gz')


# afni install dir / surfaces dir / atlas dir

source( sprintf('%s/AFNIio.R', args[2] ) )
source( sprintf('%s/coordinateFromLinearIndex.r', args[3] ) )
source( sprintf('%s/linearIndexFromCoordinate.r', args[3] ) )
library(RANN)
library(ANTsR)
mainDir <- getwd()
atlasDir <- args[4]
white_matter_mask <- args[7]
distanceThreshold <- as.numeric( strsplit( args[5], '[-]' )[[1]] ) #put here an array of 7 values
if (length(distanceThreshold)!=6) { # missing argument
  msg <- sprintf( 'distance array should be of 6 numbers' )
  warning( msg )
  stopifnot(flagDir)  
}

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
if ( length( index01 )!=0 ) {
  antLimit <- origin[ index01 ]
  posLimit <- volLimit[ index01 ]
}
if ( length( index02 )!=0 ) {
  antLimit <- volLimit[ index02 ]
  posLimit <- origin[ index02 ]
}
system('rm orient.1D origin.1D dimensions.1D voxSize.1D')

limits <- round( seq( posLimit, antLimit, length.out = 7  ), 0 )

instr <- sprintf( '3dresample -master inputAnat_mask_box.nii.gz -input %s -rmode NN -prefix white_matter_automask.nii.gz', args[7] )
system( instr )
white_matter_file <- read.AFNI( 'white_matter_automask.nii.gz' )
white_matter_volume <- white_matter_file$brk[,,,1] 

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
  
  clipVolumeLoad <- read.AFNI( sprintf( 'clip_0%1.0f.nii.gz', counter ) )
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
  
  volumeAnts <- antsImageRead( sprintf( 'clip_0%1.0f_thr.nii.gz', counter ) );
  volumeAntsN3 <- n3BiasFieldCorrection( volumeAnts, 5 )
  antsImageWrite( volumeAntsN3, 'volumeN3.nii.gz' ) 
  
  instr <- '3dcalc -a volumeN3.nii.gz -b white_matter_automask.nii.gz -expr \u0027 a*within(b,0,0) \u0027 -prefix volumeN3_mask.nii.gz'
  system( instr )
  instr <- '3dcalc -a volumeN3_mask.nii.gz -expr \u0027step(a)\u0027 -prefix maskN3.nii.gz'
  system( instr )
  instr <- '3dkmeans -f volumeN3.nii.gz -k 3 -mask maskN3.nii.gz -prefix kMeansOut03.nii.gz -r 200'
  system( instr )

  instr <- '3dcalc -a volumeN3.nii.gz -expr \u0027step(a)\u0027 -prefix maskN3_01.nii.gz'
  system( instr )
  instr <- '3dkmeans -f volumeN3.nii.gz -mask maskN3_01.nii.gz -k 3 -prefix kMeansOut03.nii.gz -r 200'
  system( instr )  
  
  instr <- sprintf('3dSeg -anat volumeN3.nii.gz -mask maskVolume.nii.gz -classes \u0027CSF ; GM ; WM\u0027 -bias_classes \u0027GM ; WM\u0027 -bias_fwhm 0.0 -mixfrac UNI -main_N %s -blur_meth BFT', args[6] ) 
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
  
  anatVol <- read.AFNI( 'Anat+orig' )
  anatVol3d <- anatVol$brk[,,,1]
  classVol <- read.AFNI( 'Classes+orig' )
  classData3d <- classVol$brk[,,,1]
  emptyVol <- array( 0, dim( classData3d ) )
  emptyVol[ white_matter_volume == 1 & classData3d > 0 ] <- 1
  emptyVol[ classData3d == 2 ] <- 1
  classVol$brk[,,,1] <- emptyVol 
  write.AFNI( filename='gmMask.nii.gz',
              brk = classVol$brk,
              origin = classVol$origin, 
              delta = classVol$delta, 
              defhead = classVol$NI_head )
  system( '3dmask_tool -input gmMask.nii.gz -dilate_input -1 +1 -prefix gmMask01.nii.gz' )
  system( '3dmask_tool -input gmMask01.nii.gz -fill_holes -prefix gmMask02.nii.gz' )
  
  imgOriginal <- 'gmMask02.nii.gz'
  innerOriginal <- read.AFNI( imgOriginal )
  innerVolumeOriginal <- innerOriginal$brk[,,,1]
  
  inside <- which( white_matter_volume == 1 & classData3d > 0 )
  outside <- which( white_matter_volume == 0 & classData3d > 0 )

  emptyVol <- array( 0, dim(innerVolumeOriginal) )
  emptyVol[inside] <- round( nnOutIn$nn.dists*-1, 2 )
  emptyVol[outside] <- round( nnOut$nn.dists, 2 )
  distanceVol <- emptyVol
  distanceThr <- as.numeric( distanceVol < distanceThreshold[counter] & distanceVol > 0 )
  innerOriginal$brk[,,,1] <- distanceThr
  write.AFNI('distanceThr03.nii.gz',
             brk=innerOriginal$brk,
             label=NULL,
             view='+orig',
             orient=innerOriginal$orient,
             origin=innerOriginal$origin,
             defhead=innerOriginal$NI_head )
  
  # emptyVol <- array( 0, dim(innerVolumeOriginal) )
  # emptyVol[inside] <- 1
  # innerOriginal$brk[,,,1] <- emptyVol
  # write.AFNI('del01.nii.gz',
  #            brk=innerOriginal$brk,
  #            label=NULL,
  #            view='+orig',
  #            orient=innerOriginal$orient,
  #            origin=innerOriginal$origin,
  #            defhead=innerOriginal$NI_head )
  # 
  coordsInside <- coordinateFromLinearIndex( inside, dim(innerVolumeOriginal) )
  coordsOutside <- coordinateFromLinearIndex( outside, dim(innerVolumeOriginal) )
  
  iPart <- t( coordsInside )
  oPart <- t( coordsOutside )
  nnOut <- nn2( iPart, oPart, k=1 )
  nnOutIn <- nn2( oPart, iPart, k=1 )
  
  emptyVol <- array( 0, dim(innerVolumeOriginal) )
  emptyVol[inside] <- round( nnOutIn$nn.dists*-1, 2 )
  emptyVol[outside] <- round( nnOut$nn.dists, 2 )
  distanceVol <- emptyVol
  distanceThr <- as.numeric( distanceVol < distanceThreshold[counter] & distanceVol > 0 )
  innerOriginal$brk[,,,1] <- distanceThr
  write.AFNI('distanceThr03.nii.gz',
             brk=innerOriginal$brk,
             label=NULL,
             view='+orig',
             orient=innerOriginal$orient,
             origin=innerOriginal$origin,
             defhead=innerOriginal$NI_head )
  
  instr <- '3dcalc -a distanceThr03.nii.gz -b gmMask02.nii.gz -expr \u0027 and(a,b) \u0027 -prefix gmMask03.nii.gz'
  system( instr )
  
  # apply clustering
  instr <- '3dclust -prefix gmMask04.nii.gz 0 100 gmMask03.nii.gz'
  system( instr )
  system( '3dmask_tool -input gmMask04.nii.gz -fill_holes -prefix gmMask05.nii.gz' )
  
  
  
}

setwd(mainDir)
dirSegmentation <- dir( pattern=sprintf( 'Segsy*') )
letterArray <- letters[1:length( dirSegmentation )]
for (k in 1:length( dirSegmentation ) ) {
  if (k==1) { instrVar <- sprintf( '-%s %s/gmMask05.nii.gz', letterArray[k], dirSegmentation[k] ) }
  if (k>1) { instrVar <- sprintf( '%s -%s %s/gmMask05.nii.gz', instrVar, letterArray[k], dirSegmentation[k] ) }
  if (k==1) { instrExpr <- sprintf( '\u0027%s', letterArray[k] ) }
  if (k>1) { instrExpr <- sprintf( '%s + %s',instrExpr, letterArray[k] ) } 
  if (k==length( dirSegmentation )) { instrExpr <- sprintf( '%s\u0027', instrExpr, letterArray[k] ) } 
}
instrCalc <- sprintf('3dcalc %s -expr %s -prefix gm_segVol.nii.gz', instrVar, instrExpr)
system( instrCalc )

# isolate WM:
instr <- '3dclust 0 100 wm_segVol.nii.gz > out.1D'
system( instr )
clustTable <- read.table( 'out.1D', comment.char = "#" )
system( 'rm out.1D' )
if ( dim(clustTable)[1]==1 ) {
  clusteringInst <- sprintf('3dclust -prefix wmMask_clust.nii.gz 0 %1.0f wm_segVol.nii.gz', clustTable[2,1] - 1 )
}
if ( dim(clustTable)[1]>1 ) {
  if (clustTable[2,1]/clustTable[1,1] >= args[6]) {
    clusteringInst <- sprintf('3dclust -prefix wmMask_clust.nii.gz 0 %1.0f wm_segVol.nii.gz', clustTable[2,1] - 1 )
  }
  if (clustTable[2,1]/clustTable[1,1] < args[6]) {
    clusteringInst <- sprintf('3dclust -prefix wmMask_clust.nii.gz 0 %1.0f wm_segVol.nii.gz', clustTable[2,1] + 1 )
  }
}
system( clusteringInst )
system( '3dmask_tool -input wmMask_clust.nii.gz -fill_holes -prefix wmMask_clust_fill.nii.gz' )
system( '3dmask_tool -input wmMask_clust_fill.nii.gz -dilate_input 1 -1 -prefix wmMask_clust_fill01.nii.gz' )
system( '3dmask_tool -input wmMask_clust_fill01.nii.gz -fill_holes -prefix wmMask_clust_fill02.nii.gz' )
system( sprintf( '3dresample -inset wmMask_clust_fill02.nii.gz -prefix wmMask_clust_fill03.nii.gz -master %s -rmode NN', args[1] ) )

# clean up
system('rm atlas*')
system('rm mask*')
system('rm wmMask_clust.nii.gz')
system('rm wm_segVol.nii.gz')
system('rm -R Segsy*')
system('rm inputAnat*')
system('rm wmMask_clust_fill.nii.gz')
system('rm wmMask_clust_fill01.nii.gz')
system('rm wmMask_clust_fill02.nii.gz')
system('mv wmMask_clust_fill03.nii.gz white_matter_mask.nii.gz')
