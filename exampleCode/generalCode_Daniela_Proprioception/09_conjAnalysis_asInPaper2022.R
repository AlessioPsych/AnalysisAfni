rm(list=ls())
atlasFolder <- '/analyse/Project0226/tests/Daniela_pilotData/suma_MNI152_2009_princetonAtlas'
setwd( atlasFolder )
getwd()

thrValue <- 4.5
thrValueMask <- 2.7


#### linear non REML ####

if ( file.exists( sprintf( 'tOutcome_mask_linear_asInPaper2022.nii.gz' ) ) ) {
  instr <- sprintf( 'rm tOutcome_mask_linear_asInPaper2022.nii.gz' ); system( instr )
} 

if ( file.exists( sprintf( 'tOutcome_linear_conjunction_asInPaper2022.nii.gz' ) ) ) {
  instr <- sprintf( 'rm tOutcome_linear_conjunction_asInPaper2022.nii.gz' ); system( instr )
} 

instr <- paste( '3dcalc',
                '-a tOutcome_linear_23.nii.gz',
                sprintf( '-expr \u027 not( step(a-%1.4f) ) \u027', thrValueMask ),
                '-prefix tOutcome_mask_linear_asInPaper2022.nii.gz' )
print( instr )
system( instr )

instr <- paste( '3dcalc',
                '-a tOutcome_linear_17.nii.gz',
                '-b tOutcome_linear_20.nii.gz',
                '-c tOutcome_linear_mask_asInPaper2022.nii.gz',
                sprintf( '-expr \u027 step(a-%1.4f)*c + step(b-%1.4f)*c \u027', thrValue, thrValue ),
                '-prefix tOutcome_linear_conjunction_asInPaper2022.nii.gz' )
print( instr )
system( instr )

#### non-linear non REML ####

if ( file.exists( sprintf( 'tOutcome_mask_nonlinear_asInPaper2022.nii.gz' ) ) ) {
  instr <- sprintf( 'rm tOutcome_mask_nonlinear_asInPaper2022.nii.gz' ); system( instr )
} 

if ( file.exists( sprintf( 'tOutcome_nonlinear_conjunction_asInPaper2022.nii.gz' ) ) ) {
  instr <- sprintf( 'rm tOutcome_nonlinear_conjunction_asInPaper2022.nii.gz' ); system( instr )
} 

instr <- paste( '3dcalc',
                '-a tOutcome_nonlinear_23.nii.gz',
                sprintf( '-expr \u027 not( step(a-%1.4f) ) \u027', thrValueMask ),
                '-prefix tOutcome_mask_nonlinear_asInPaper2022.nii.gz' )
print( instr )
system( instr )

instr <- paste( '3dcalc',
                '-a tOutcome_nonlinear_17.nii.gz',
                '-b tOutcome_nonlinear_20.nii.gz',
                '-c tOutcome_nonlinear_mask_asInPaper2022.nii.gz',
                sprintf( '-expr \u027 step(a-%1.4f)*c + step(b-%1.4f)*c \u027', thrValue, thrValue ),
                '-prefix tOutcome_nonlinear_conjunction_asInPaper2022.nii.gz' )
print( instr )
system( instr )

#### linear REML ####

if ( file.exists( sprintf( 'tOutcome_linear_REML_mask_asInPaper2022.nii.gz' ) ) ) {
  instr <- sprintf( 'rm tOutcome_linear_REML_mask_asInPaper2022.nii.gz' ); system( instr )
} 

if ( file.exists( sprintf( 'tOutcome_linear_REML_conjunction_asInPaper2022.nii.gz' ) ) ) {
  instr <- sprintf( 'rm tOutcome_linear_REML_conjunction_asInPaper2022.nii.gz' ); system( instr )
} 

instr <- paste( '3dcalc',
                '-a tOutcome_REML_linear_23.nii.gz',
                sprintf( '-expr \u027 not( step(a-%1.4f) ) \u027', thrValueMask ),
                '-prefix tOutcome_linear_REML_mask_asInPaper2022.nii.gz' )
print( instr )
system( instr )

instr <- paste( '3dcalc',
                '-a tOutcome_REML_linear_17.nii.gz',
                '-b tOutcome_REML_linear_20.nii.gz',
                '-c tOutcome_linear_REML_mask_asInPaper2022.nii.gz',
                sprintf( '-expr \u027 step(a-%1.4f)*c + step(b-%1.4f)*c \u027', thrValue, thrValue ),
                '-prefix tOutcome_linear_REML_conjunction_asInPaper2022.nii.gz' )
print( instr )
system( instr )

#### non linear REML ####

if ( file.exists( sprintf( 'tOutcome_nonlinear_REML_mask_asInPaper2022.nii.gz' ) ) ) {
  instr <- sprintf( 'rm tOutcome_nonlinear_REML_mask_asInPaper2022.nii.gz' ); system( instr )
} 

if ( file.exists( sprintf( 'tOutcome_nonlinear_REML_conjunction_asInPaper2022.nii.gz' ) ) ) {
  instr <- sprintf( 'rm tOutcome_nonlinear_REML_conjunction_asInPaper2022.nii.gz' ); system( instr )
} 

instr <- paste( '3dcalc',
                '-a tOutcome_REML_nonlinear_23.nii.gz',
                sprintf( '-expr \u027 not( step(a-%1.4f) ) \u027', thrValueMask ),
                '-prefix tOutcome_nonlinear_REML_mask_asInPaper2022.nii.gz' )
print( instr )
system( instr )

instr <- paste( '3dcalc',
                '-a tOutcome_REML_nonlinear_17.nii.gz',
                '-b tOutcome_REML_nonlinear_20.nii.gz',
                '-c tOutcome_nonlinear_REML_mask_asInPaper2022.nii.gz',
                sprintf( '-expr \u027 step(a-%1.4f)*c + step(b-%1.4f)*c \u027', thrValue, thrValue ),
                '-prefix tOutcome_nonlinear_REML_conjunction_asInPaper2022.nii.gz' )
print( instr )
system( instr )



