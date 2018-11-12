args <- commandArgs(T)
print( args )

#setwd('/media/alessiofracasso/storage2/dataLaminar_leftVSrightNEW/V2677leftright_Ric/resultsCBS')
#args <- c( '-node_normals' )

surfaceFilesTopo <- dir( 'surfaces_folder/', pattern='*.topo'  )
surfaceFilesTopo <- surfaceFilesTopo[1:length(surfaceFilesTopo)-1]
surfaceFilesCoords <- dir( 'surfaces_folder/', pattern='*.coord'  )
surfaceFilesCoords <- surfaceFilesCoords[1:length(surfaceFilesCoords)-1]
surfaceFileSpec <- dir( 'surfaces_folder/', pattern='spec*'  )

dirSplit <- strsplit(args[1],'[-]')
newDir <- sprintf('surface_%s', dirSplit[[1]][2] )

## checkdir, create if it does not exists, if it exists then halts execution
mainDir <- getwd()
flagDir <- dir.create( file.path(mainDir, newDir) )
if (flagDir==FALSE) { # directory already exists!
  msg <- sprintf( 'Remove the directory %s_folder to proceed', args[1] )
  warning( msg )
  stopifnot(flagDir)  
}

for ( k in 1:length(surfaceFilesCoords) ) {
  outSplitName <- strsplit(surfaceFilesCoords[k],'[.]')
  outputFile <- sprintf( '%s/%s.normals', newDir, outSplitName[[1]][1] )
  instr <- sprintf('SurfaceMetrics %s -i_1D surfaces_folder/%s surfaces_folder/%s -prefix %s',
                   args[1], surfaceFilesCoords[k], surfaceFilesTopo[k], outputFile )  
  if (k<10) { system( sprintf( 'echo surface n. 0%s, files: %s %s', k, surfaceFilesCoords[k], surfaceFilesTopo[k] ) ) }
  if (k>=10) { system( sprintf( 'echo surface n. %s, files: %s %s', k, surfaceFilesCoords[k], surfaceFilesTopo[k] ) ) }
  
  system( instr )
  
}

#SurfSmooth  -spec quick.spec -surf_A NodeList.1D -met HEAT_05   \
#-input in.1D -fwhm 8 -add_index         \
#-output in_smh8.1D.dset 