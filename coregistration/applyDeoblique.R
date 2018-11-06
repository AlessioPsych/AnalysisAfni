args <- commandArgs(T)
print( args )

#args <- c('niftiFiles_delme/','nifti_deob/')

inputDir <- args[1]
outDir <- args[2]
mainDir <- getwd()
setwd( mainDir )

targetDir <- outDir
flagDir <- dir.create( file.path(mainDir, targetDir) )
if (flagDir==FALSE) { # directory already exists!
  msg <- sprintf( 'Remove the directory %s to proceed', targetDir )
  warning( msg )
  stopifnot(flagDir)  
}


allFiles <- list.files( inputDir )
filesNii <- allFiles[ grep('.nii',allFiles,fixed=T) ]

setwd( mainDir )
for (k in 1:length(filesNii)) {
  if (k < 10) { instr <- sprintf('3dWarp -deoblique -prefix %s0%s_deob.nii %s%s',
                                 inputDir, k, inputDir, filesNii[k] ) }
  if (k >= 10) { instr <- sprintf('3dWarp -deoblique -prefix %s%s_deob.nii %s%s',
                                  inputDir, k, inputDir, filesNii[k] ) }
  print( instr )
  system( instr )
}

system( sprintf('mv %s*deob* %s', inputDir, targetDir) )
setwd( mainDir )

#system('rm _ttt*')
#system('rm *.nii')
#system('gunzip *.gz')
