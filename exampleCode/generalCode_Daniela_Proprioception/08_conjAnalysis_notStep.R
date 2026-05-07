rm(list=ls())
atlasFolder <- '/analyse/Project0226/tests/Daniela_pilotData/suma_MNI152_2009_princetonAtlas'
setwd( atlasFolder )
getwd()

if ( file.exists( sprintf( 'tOutcome_linear_conjunction_notStep.nii.gz' ) ) ) {
  instr <- sprintf( 'rm tOutcome_linear_conjunction_notStep.nii.gz' ); system( instr )
} 
instr <- paste( '3dcalc',
                '-prefix conjunction_map',
                '-a tOutcome_linear_17.nii.gz',
                '-b tOutcome_linear_20.nii.gz', 
                '-c tOutcome_linear_23.nii.gz', 
                '-d tOutcome_linear_26.nii.gz', 
                sprintf( '-expr \u027 a+b-c-d  \u027' ),
                '-prefix tOutcome_linear_conjunction_notStep.nii.gz' )
print( instr )
system( instr )

if ( file.exists( sprintf( 'tOutcome_nonlinear_conjunction_notStep.nii.gz' ) ) ) {
  instr <- sprintf( 'rm tOutcome_nonlinear_conjunction_notStep.nii.gz' ); system( instr )
} 
instr <- paste( '3dcalc',
                '-prefix conjunction_map',
                '-a tOutcome_nonlinear_17.nii.gz',
                '-b tOutcome_nonlinear_20.nii.gz', 
                '-c tOutcome_nonlinear_23.nii.gz', 
                '-d tOutcome_nonlinear_26.nii.gz', 
                sprintf( '-expr \u027 a+b-c-d  \u027' ),
                '-prefix tOutcome_nonlinear_conjunction_notStep.nii.gz' )
print( instr )
system( instr )

if ( file.exists( sprintf( 'tOutcome_REML_linear_conjunction_notStep.nii.gz' ) ) ) {
  instr <- sprintf( 'rm tOutcome_REML_linear_conjunction_notStep.nii.gz' ); system( instr )
} 
instr <- paste( '3dcalc',
                '-prefix conjunction_map',
                '-a tOutcome_REML_linear_17.nii.gz',
                '-b tOutcome_REML_linear_20.nii.gz', 
                '-c tOutcome_REML_linear_23.nii.gz', 
                '-d tOutcome_REML_linear_26.nii.gz', 
                sprintf( '-expr \u027 a+b-c-d  \u027' ),
                '-prefix tOutcome_REML_linear_conjunction_notStep.nii.gz' )
print( instr )
system( instr )

if ( file.exists( sprintf( 'tOutcome_REML_nonlinear_conjunction_notStep.nii.gz' ) ) ) {
  instr <- sprintf( 'rm tOutcome_REML_nonlinear_conjunction_notStep.nii.gz' ); system( instr )
} 
instr <- paste( '3dcalc',
                '-prefix conjunction_map',
                '-a tOutcome_REML_nonlinear_17.nii.gz',
                '-b tOutcome_REML_nonlinear_20.nii.gz', 
                '-c tOutcome_REML_nonlinear_23.nii.gz', 
                '-d tOutcome_REML_nonlinear_26.nii.gz', 
                sprintf( '-expr \u027 a+b-c-d  \u027' ),
                '-prefix tOutcome_REML_nonlinear_conjunction_notStep.nii.gz' )
print( instr )
system( instr )
