args <- commandArgs(T)
print( args )

#instr01 <- sprintf( '3dWarp -card2oblique %s -prefix ts.nii  %s', args[1], args[2] ) 

instr02 <- sprintf( '@Align_Centers -cm -dset %s -base %s', args[1], args[2] )

arg1split <- strsplit(args[1],'.nii')
shiftFileName1D <- sprintf( '%s_shft.1D', arg1split[[1]][1] )
shiftFileName1D
shiftFileName <- sprintf( '%s_shft.nii', arg1split[[1]][1] )
shiftFileName

outputFileName <- sprintf( '%s_t.nii', arg1split[[1]][1] )
outputFileName

#argOut <- args[1]
#length(as.character(argOut))
#argOut[1:length(argOut)-3]
#outputFileName

instr03 <- sprintf( 'mri_convert --apply_transform %s --in_type nii --out_type nii %s %s', args[3], shiftFileName,  outputFileName)

#instr01
instr02
instr03

#system( instr01 )
system( instr02 )
system( instr03 )

system( sprintf( 'rm %s',  shiftFileName1D ) )
system( sprintf( 'rm %s',  shiftFileName ) )
