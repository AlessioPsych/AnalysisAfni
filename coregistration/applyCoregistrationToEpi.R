args <- commandArgs(T)
print( args )

#args <- c('nifti_deob/','coregistration/meanEpiSingleShot.nii.gz','coregistration/tMat.1D','nifti_coreg/')

dirFiles <- args[1]
coregTarget <- args[2]
coregMat <- args[3]
outDir <- args[4]
mainDir <- getwd()
setwd( mainDir )

targetDir <- outDir
flagDir <- dir.create( file.path(mainDir, targetDir) )
if (flagDir==FALSE) { # directory already exists!
  msg <- sprintf( 'Remove the directory %s to proceed', targetDir )
  warning( msg )
  stopifnot(flagDir)  
}


allFiles <- list.files( dirFiles )
filesNii <- allFiles[ grep('.nii',allFiles,fixed=T) ]

instr <- sprintf( '3dAutobox -noclust -prefix %s_ttt_TargetVol.nii %s', dirFiles, coregTarget )
system( instr )

setwd( mainDir )
for (k in 1:length(filesNii)) {
  if (k < 10) { instr <- sprintf('3dAllineate -final linear -master %s_ttt_TargetVol.nii -source %s%s -1Dmatrix_apply %s -prefix %s0%s_coreg.nii',
                                 dirFiles, dirFiles,filesNii[k], coregMat, dirFiles, k ) }
  if (k >= 10) { instr <- sprintf('3dAllineate -final linear -master %s_ttt_TargetVol.nii -source %s%s -1Dmatrix_apply %s -prefix %s%s_coreg.nii',
                                  dirFiles, dirFiles,filesNii[k], coregMat, dirFiles, k ) }
  print( instr )
  system( instr )
}

instr <- sprintf( 'rm %s_ttt_TargetVol.nii', dirFiles )
system( instr )

system( sprintf('mv %s*coreg* %s', dirFiles, targetDir) )

#system('rm _ttt*')
#system('rm *.nii')
#system('gunzip *.gz')
#setwd( mainDir )
