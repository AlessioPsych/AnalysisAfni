rm(list=ls())
atlasFolder <- '/analyse/Project0226/tests/Daniela_pilotData/suma_MNI152_2009_princetonAtlas'
setwd( atlasFolder )
getwd()

thrValue <- 2.7

if ( file.exists( sprintf( 'tOutcome_linear_conjunction.nii.gz' ) ) ) {
  instr <- sprintf( 'rm tOutcome_linear_conjunction.nii.gz' ); system( instr )
} 
instr <- paste( '3dcalc',
                '-prefix conjunction_map',
                '-a tOutcome_linear_17.nii.gz',
                '-b tOutcome_linear_20.nii.gz', 
                '-c tOutcome_linear_23.nii.gz', 
                '-d tOutcome_linear_26.nii.gz', 
                sprintf( '-expr \u027 step(a-%1.4f) + step(b-%1.4f) - step(c-%1.4f) - step(d-%1.4f)  \u027', thrValue, thrValue, thrValue, thrValue ),
                '-prefix tOutcome_linear_conjunction.nii.gz' )
print( instr )
system( instr )

if ( file.exists( sprintf( 'tOutcome_nonlinear_conjunction.nii.gz' ) ) ) {
  instr <- sprintf( 'rm tOutcome_nonlinear_conjunction.nii.gz' ); system( instr )
} 
instr <- paste( '3dcalc',
                '-prefix conjunction_map',
                '-a tOutcome_nonlinear_17.nii.gz',
                '-b tOutcome_nonlinear_20.nii.gz', 
                '-c tOutcome_nonlinear_23.nii.gz', 
                '-d tOutcome_nonlinear_26.nii.gz', 
                sprintf( '-expr \u027 step(a-%1.4f) + step(b-%1.4f) - step(c-%1.4f) - step(d-%1.4f)  \u027', thrValue, thrValue, thrValue, thrValue ),
                '-prefix tOutcome_nonlinear_conjunction.nii.gz' )
print( instr )
system( instr )

if ( file.exists( sprintf( 'tOutcome_REML_linear_conjunction.nii.gz' ) ) ) {
  instr <- sprintf( 'rm tOutcome_REML_linear_conjunction.nii.gz' ); system( instr )
} 
instr <- paste( '3dcalc',
                '-prefix conjunction_map',
                '-a tOutcome_REML_linear_17.nii.gz',
                '-b tOutcome_REML_linear_20.nii.gz', 
                '-c tOutcome_REML_linear_23.nii.gz', 
                '-d tOutcome_REML_linear_26.nii.gz', 
                sprintf( '-expr \u027 step(a-%1.4f) + step(b-%1.4f) - step(c-%1.4f) - step(d-%1.4f)  \u027', thrValue, thrValue, thrValue, thrValue ),
                '-prefix tOutcome_REML_linear_conjunction.nii.gz' )
print( instr )
system( instr )

if ( file.exists( sprintf( 'tOutcome_REML_nonlinear_conjunction.nii.gz' ) ) ) {
  instr <- sprintf( 'rm tOutcome_REML_nonlinear_conjunction.nii.gz' ); system( instr )
} 
instr <- paste( '3dcalc',
                '-prefix conjunction_map',
                '-a tOutcome_REML_nonlinear_17.nii.gz',
                '-b tOutcome_REML_nonlinear_20.nii.gz', 
                '-c tOutcome_REML_nonlinear_23.nii.gz', 
                '-d tOutcome_REML_nonlinear_26.nii.gz', 
                sprintf( '-expr \u027 step(a-%1.4f) + step(b-%1.4f) - step(c-%1.4f) - step(d-%1.4f)  \u027', thrValue, thrValue, thrValue, thrValue ),
                '-prefix tOutcome_REML_nonlinear_conjunction.nii.gz' )
print( instr )
system( instr )
