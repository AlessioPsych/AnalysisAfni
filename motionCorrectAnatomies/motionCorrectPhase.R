args <- commandArgs(T)
print( args )

#setwd('/media/alessiofracasso/storage2/highResMyelinT1/wholeBrainHighResMyelinAnatomies/Jeroen_myelinstuff/afni')
#args <- c( 'phase/', '*crop.nii.gz', '_al_mat.aff12.1D' )

mainDir <- getwd()
selDir <- args[1]

# get in the directory and select filenames (volumes + matrices)
setwd( sprintf('%s/%s', mainDir, selDir ) )
filesVolume <- dir( pattern=args[2] )
filesMat <- dir( pattern=args[3] )

# apply motion correction 
for ( k in 1 : length(filesVolume) ) {
#for ( k in 1 : 1 ) {
  
  splitFileName <- strsplit( filesVolume[k], '[.]' )
  splitFileNameOut <- sprintf( '%s_al.nii.gz', splitFileName[[1]][1] )
  
  if (k == 1) {
    
    write.table( t( rep(0,12) ), file='tempParamFile.1D', col.names = FALSE, row.names = FALSE )
    
    instr01 <- sprintf( '3dAllineate -source %s -master %s -prefix %s -final wsinc5 -1Dparam_apply tempParamFile.1D', filesVolume[k] ,filesVolume[k], splitFileNameOut )
    
    system( sprintf('echo %s', instr01) )
    
    system( instr01 )
    
    system( 'rm tempParamFile.1D' )
    
  }
  
  if (k > 1) {
    
    instr01 <- sprintf( '3dAllineate -1Dmatrix_apply %s -source %s -prefix %s -master %s -final wsinc5', filesMat[k-1], filesVolume[k], splitFileNameOut, filesVolume[k] )
    
    system( sprintf('echo %s', instr01) )
    
    system( instr01 )
    
  }
  
}

setwd( sprintf('%s', mainDir ) )