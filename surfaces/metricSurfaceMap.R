args <- commandArgs(T)
print( args )

inputMetric <- args[1]
inputSurfaceFolder <- args[2]

#setwd('/media/alessiofracasso/storage2/dataLaminar_leftVSrightNEW/V2677leftright_Ric/resultsCBS')
#args <- c( '-node_normals' )

surfaceFilesTopo <- dir( inputSurfaceFolder, pattern='*.topo'  )
surfaceFilesTopo <- surfaceFilesTopo[1:length(surfaceFilesTopo)-1]
surfaceFilesCoords <- dir( inputSurfaceFolder, pattern='*.coord'  )
surfaceFilesCoords <- surfaceFilesCoords[1:length(surfaceFilesCoords)-1]
surfaceFileSpec <- dir( inputSurfaceFolder, pattern='spec*'  )

dirSplit <- strsplit( inputMetric,'[-]' )
newDir <- sprintf('%s_%s', inputSurfaceFolder, dirSplit[[1]][2] )

## checkdir, create if it does not exists, if it exists then halts execution
mainDir <- getwd()
flagDir <- dir.create( file.path(mainDir, newDir) )
if (flagDir==FALSE) { # directory already exists!
  msg <- sprintf( 'Remove the directory %s to proceed', newDir )
  warning( msg )
  stopifnot(flagDir)  
}

for ( k in 1:length(surfaceFilesCoords) ) {
  outSplitName <- strsplit(surfaceFilesCoords[k],'[.]')
  outputFile <- sprintf( '%s/%s.metric', newDir, outSplitName[[1]][1] )
  instr <- sprintf('SurfaceMetrics %s -i_1D %s/%s %s/%s -prefix %s',
                   inputMetric, inputSurfaceFolder, surfaceFilesCoords[k], inputSurfaceFolder, surfaceFilesTopo[k], outputFile )  
  if (k<10) { system( sprintf( 'echo surface n. 0%s, files: %s %s', k, surfaceFilesCoords[k], surfaceFilesTopo[k] ) ) }
  if (k>=10) { system( sprintf( 'echo surface n. %s, files: %s %s', k, surfaceFilesCoords[k], surfaceFilesTopo[k] ) ) }
  
  system( instr )
  
}

#SurfSmooth  -spec quick.spec -surf_A NodeList.1D -met HEAT_05   \
#-input in.1D -fwhm 8 -add_index         \
#-output in_smh8.1D.dset 