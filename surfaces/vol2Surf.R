args <- commandArgs(T)
print( args )

mainDir <- getwd()


#rm(list=ls())
#mainDir <- '/home/fracasso/data/ArjanProject/Arjan_fMRI_scripts/results'
#setwd( mainDir )
#args <- c('leftSmall.1D.roi_clust+orig')


roiName <- strsplit( args[1], '[.]'  )[[1]][1]

instr <- sprintf('3dcopy %s ttt_roiVolume.nii.gz', args[1])
system( instr )

setwd('surfaces_folder')
boundaries <- dir(pattern='.*boundary.*coord')

setwd( mainDir )
for (bCounter in 1:length(boundaries) ) {
  print( sprintf('boundary: %s, file: %s', as.character(bCounter), boundaries[bCounter] ) )
  surfaceName <- boundaries[bCounter]
  surfOutputName <- ifelse( bCounter<10, 
                      sprintf('roiProjOnSurface_boundary0%s.1D.dset', as.character(bCounter-1) ), 
                      sprintf('roiProjOnSurface_boundary%s.1D.dset', as.character(bCounter-1) ) )
  instr <- c( '3dVol2Surf', '-spec surfaces_folder/spec.surfaces.smoothed',
               sprintf('-surf_A surfaces_folder/%s', surfaceName),
               '-sv ttt_roiVolume.nii.gz',
               '-grid_parent ttt_roiVolume.nii.gz',
               '-map_func mask',
               sprintf( '-out_1D %s', surfOutputName ) )
  instrConcatenated <- paste( instr, collapse = " "  )
  system( instrConcatenated )
}

setwd( mainDir )
subDir <- sprintf( '%s_surfaceDir', roiName )
flagDir <- dir.create( file.path(mainDir, subDir) )
if (flagDir==FALSE) { # directory already exists!
  msg <- sprintf( 'Remove the directory %s to proceed', subDir )
  warning( msg )
  stopifnot(flagDir)  
}

mvInstr <- sprintf( 'mv roiProjOnSurface_boundary* %s_surfaceDir', roiName )
system( mvInstr )

system('rm ttt_roiVolume.nii.gz')

