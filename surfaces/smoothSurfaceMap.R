args <- commandArgs(T)
print( args )

#setwd('/media/alessiofracasso/storage2/dataLaminar_leftVSrightNEW/V2677leftright_Ric/resultsCBS')
#args <- c( 'curvatureCopy_interp_folder', '8' )
#print( args )

surfaceFilesTopo <- dir( 'surfaces_folder/', pattern='*.topo'  )
surfaceFilesTopo <- surfaceFilesTopo[1:length(surfaceFilesTopo)-1]
surfaceFilesCoords <- dir( 'surfaces_folder/', pattern='*.coord'  )
surfaceFilesCoords <- surfaceFilesCoords[1:length(surfaceFilesCoords)-1]
surfaceFileSpec <- dir( 'surfaces_folder/', pattern='spec*'  )
mapfiles <- dir( args[1] )

dirSplit <- strsplit(args[1],'[/]')
newDir <- sprintf('%s_smoothed_%s', dirSplit[[1]][1], args[2] )

## checkdir, create if it does not exists, if it exists then halts execution
mainDir <- getwd()
flagDir <- dir.create( file.path(mainDir, newDir) )
if (flagDir==FALSE) { # directory already exists!
  msg <- sprintf( 'Remove the directory %s_folder to proceed', args[1] )
  warning( msg )
  stopifnot(flagDir)  
}

for ( k in 1:length(surfaceFilesCoords) ) {
  outSplitName <- strsplit(mapfiles[k],'[.]')
  newDirSplit <- strsplit(newDir,'[_]')
  outputFile <- sprintf( '%s/%s_%s_%s_smooth', newDir, outSplitName[[1]][1], newDirSplit[[1]][1], args[2] )
  instr <- sprintf('SurfSmooth -spec surfaces_folder/%s -surf_A surfaces_folder/%s -met HEAT_05 -input %s%s -fwhm %s -output %s',
                   surfaceFileSpec, surfaceFilesCoords[k], args[1], mapfiles[k], args[2], outputFile )  
  system( instr )
}

#SurfSmooth  -spec quick.spec -surf_A NodeList.1D -met HEAT_05   \
#-input in.1D -fwhm 8 -add_index         \
#-output in_smh8.1D.dset 
