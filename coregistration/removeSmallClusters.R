args <- commandArgs(T)
print(args)
#args <- c('_ttt_T1MPRAGE_SS_N4.nii.gz', '/home/alessiof/abin', '0')
#args <- c('pdCorrectRegularT1_noBlur_stripped.nii','/packages/afni/17.0.13','1')
#'TT_152_2009c+tlrc'
#setwd('/home/dijkj/Documents/Segmentation_tests/JeDi/RemoveCerebellum')
#setwd('/home/alessiof/Mount/Data/anatData/anatAnalysis_06032018/20171101_VBE10_preproc')

inputAnat <- args[1]
minProportion <- as.numeric( args[2] )
outFile <- args[3]
print( inputAnat )
print( minProportion )
print( outFile )

instr <- sprintf( '3dcopy %s ttt_inputData.nii.gz', inputAnat )

system( instr )
instr <- sprintf('3dclust 0 5 ttt_inputData.nii.gz > ttt_out.1D' )
system( instr )
dataIn <- read.table( 'ttt_out.1D', as.is=TRUE )
pTable <- round( dataIn[,1] / sum( dataIn[,1]), 3 )
pIdx <- pTable<minProportion
cLimit <- dataIn[which( pIdx )[1],1]
instr <- sprintf('3dclust -prefix %s 0 %1.0f ttt_inputData.nii.gz', outFile, cLimit+1 )
system( instr )


system('rm ttt_inputData.nii.gz')
system('rm ttt_out.1D')
