args <- commandArgs(T)
print( args )

surfacesDir <- args[1]
anatomy <- args[2]
instr <- sprintf( 'afni -niml & suma -spec %s/spec.surfaces.smoothed -sv %s &', surfacesDir, anatomy )
print( instr )
system( instr )