args <- commandArgs(T)
print( args )

source( sprintf('%s/AFNIio.R', args[4] ) )
source( sprintf('%s/coordinateFromLinearIndex.r', args[5] ) )
source( sprintf('%s/linearIndexFromCoordinate.r', args[5] ) )
library(RANN)

#ROINAME='leftV1.1D.roi'
#BOUNDARYWM='boundary00'
#BOUNDARYGM='boundary01'
#ANAT='MPRAGE_al_epi.nii.gz'
#DEPTH='DEPTH_al_epi.nii.gz'
#args <- c( ROINAME, BOUNDARYWM, DEPTH, '/usr/lib/afni/bin/','/home/alessiofracasso/Dropbox/analysisAfni/surfaces/')

#get voxels on WM border
#instr <- sprintf( 'cp surfaces_folder/%s_sm.1D.coord surfaceCoords.1D', args[2] ) 
#print( instr )
#system( instr )
minClustSize <- args[6]

commandLine <- sprintf( '3dSurf2Vol -spec %sspec.surfaces.smoothed -surf_A %s%s_sm.1D.coord -sv %s -grid_parent %s -map_func ave -prefix %s -sdata_1D %s', args[7], args[7], args[2], args[3], args[3], 'surfVolWMRoi', args[1] )


print( commandLine )
system( commandLine )

system( 'rm surfaceCoords.1D' )

wmBorderMask <- read.AFNI('surfVolWMRoi+orig')
gmLevelMask <- read.AFNI( args[3] )

maskData <- wmBorderMask$brk
gmData <- gmLevelMask$brk[,,,1]

emptyVol <- maskData
stepLimit <- 0.1
limit1 <- seq(0,1-stepLimit,by=stepLimit)
limit2 <- limit1 + stepLimit
for (nVols in 1:dim(emptyVol)[4]) {
  emptyVolTemp <- emptyVol[,,,nVols]
  for (lim in 1:length(limit1) ) {
    ind1 <- which( abs(emptyVolTemp)>0.000001 )
    ind2 <- which( gmData>limit1[lim] & gmData<=limit2[lim] )
    coordsWM <- matrix( t( coordinateFromLinearIndex( ind1, dim(emptyVolTemp) ) ), ncol=3 )
    coordsGM <- matrix( t( coordinateFromLinearIndex( ind2, dim(emptyVolTemp) ) ), ncol=3 )
    nnOut <- nn2( coordsGM, coordsWM,  k=1 ) #distance from ind2 to ind1
    coords <- coordsGM[ nnOut$nn.idx, ] #coords in ind2 close to ind1
    coordsVol <- linearIndexFromCoordinate( t( coords ), dim( emptyVolTemp ) ) #index in ind2 close to ind1
    emptyVolTemp[ coordsVol ] <- emptyVolTemp[ ind1 ]
    #emptyVol[ coordsVol ] <- emptyVol[ ind1 ]
    #for (l in 1:dim(coords)[1]) {
    #  emptyVol[ coords[l,1], coords[l,2], coords[l,3] ] <- 1
    #}
  }
  emptyVol[,,,nVols] <- emptyVolTemp
}
#image(emptyVol[,,14])

roiFileName <- sprintf( '%s_boundaryOnlyFile+orig',args[1] )
write.AFNI(roiFileName,
           brk=emptyVol,
           label=NULL,
           view='+orig',
           orient=wmBorderMask$orient,
           origin=wmBorderMask$origin,
           delta=wmBorderMask$delta,
           defhead=wmBorderMask$NI_head )

#roiFileNameClust <- sprintf( '%s_clust+orig',args[1] )
clusteringInst <- sprintf('3dclust -prefix %s_clust.nii.gz 0 %s %s', args[1], minClustSize, roiFileName)
print( clusteringInst )
system( clusteringInst )

#instr <- sprintf( '3dmask_tool -input _ttt_clust.nii.gz -prefix %s -fill_holes', roiFileNameClust, roiFileName )
#print( instr )
#system( instr )

#system( sprintf('rm %s.BRIK.gz', roiFileName ) )
#system( sprintf('rm %s.HEAD', roiFileName ) )
#system('rm _ttt_clust.nii.gz')
#system('rm surfVolWMRoi+orig.HEAD')
#system('rm surfVolWMRoi+orig.BRIK.gz')
