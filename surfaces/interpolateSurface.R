args <- commandArgs(T)
print( args )

## debug:
#setwd('/media/alessiofracasso/storage2/prfSizeDepth_2/V3814Serge/anatomy_CBS')
#args<-c('anatomy_crop.nii.gz')

argSplit <- strsplit(args[1],'[.]')
fileBase <- argSplit[[1]][1]
surfacesDir <- args[2]

## checkdir, create if it does not exists, if it exists then halts execution
mainDir <- getwd()
subDir <- sprintf( '%s_interp_%s', fileBase, surfacesDir )
flagDir <- dir.create( file.path(mainDir, subDir) )
if (flagDir==FALSE) { # directory already exists!
  msg <- sprintf( 'Remove the directory %s_folder to proceed', args[1] )
  warning( msg )
  stopifnot(flagDir)  
}

specFile <- sprintf('%s/spec.surfaces.smoothed',surfacesDir)
setwd( sprintf( '%s/%s', mainDir, surfacesDir)  )
surfaceFiles <- dir(pattern='^boundary.*_sm.1D.coord')
setwd( mainDir )

for (k in 1:length( surfaceFiles) ) {
    
  surfFileSplit <- strsplit(surfaceFiles[k],'[.]')
  
  outputFileName <- sprintf('%s_%s%s', surfFileSplit[[1]][1], fileBase, '_surf.1D.dset' )
  
  instr <- sprintf('3dVol2Surf -spec %s -surf_A %s/%s -sv %s -grid_parent %s -map_func mask -out_1D %s', 
                   specFile, surfacesDir, surfaceFiles[k], args[1], args[1], outputFileName)
  
  system(instr)
  print( instr )
  print( surfaceFiles[k] )
}

## move all the files generated into the subDir
system( sprintf('mv boundary*_surf.1D.dset %s', subDir) )
