args <- commandArgs(T)
print( args )

#setwd('/media/alessiofracasso/storage2/highResMyelinT1/wholeBrainHighResMyelinAnatomies/Jeroen_myelinstuff/afni')
#args <- c( 'phase/', '*t001-0001.nii', '-anterior_25_-left_-5', '_crop' )

mainDir <- getwd()
selDir <- args[1]

setwd( sprintf('%s/%s', mainDir, selDir ) )
files <- dir( pattern=args[2] )
clipInstr <- strsplit( args[3], '_' )
for ( k in 1 : length( clipInstr[[1]]  )  ) {
  if (k==1) { strOut <- clipInstr[[1]][k] }
  if (k>1) { strOut <- sprintf( '%s %s', strOut, clipInstr[[1]][k] )   }
}

for ( k in 1 : length(files) ) {
  
  splitFileName <- strsplit( files[k], '[.]' )
  splitFileNameOut <- sprintf( '%s_crop.nii.gz', splitFileName[[1]][1] )
  
  if (k == 1) {
    
    instr01 <- sprintf( '@clip_volume %s -prefix clipped.nii.gz -input %s', strOut, files[k] )
    instr02 <- sprintf( '3dAutobox -input clipped.nii.gz -prefix clippedBox.nii.gz' )
    instr03 <- sprintf( '3dcopy clippedBox.nii.gz %s', splitFileNameOut )
    
    system( sprintf('echo %s', instr01) )
    system( sprintf('echo %s', instr02) )
    system( sprintf('echo %s', instr03) )
    
    system( instr01 )
    system( instr02 )
    system( instr03 )
    
  }
  
  if (k > 1) {
    instr01 <- sprintf( '3dresample -master clippedBox.nii.gz -inset %s -prefix %s -rmode NN', files[k], splitFileNameOut )
    system( sprintf('echo %s', instr01) )
    system( instr01 )
    
  }

}

system('rm clipped.nii.gz')
system('rm clippedBox.nii.gz')

setwd( sprintf('%s', mainDir ) )