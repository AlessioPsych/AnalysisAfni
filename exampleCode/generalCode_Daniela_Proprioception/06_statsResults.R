rm(list=ls())
atlasFolder <- '/analyse/Project0226/tests/Daniela_pilotData/suma_MNI152_2009_princetonAtlas'
setwd( atlasFolder )
getwd()

nBrik <- 26

# 2. t stat active
# 5. t stat passive
# 8. t stat touch
# 14. t stat instruction

# 17. t stat active vs rest
# 20. t stat passive vs rest
# 23. t stat touch vs rest
# 26. t stat instruction vs rest

for ( nBrik in c(2,5,8,14,17,20,23,26 ) ) {

# clean up
if ( file.exists( sprintf( 'tTestMask.nii.gz' ) ) ) {
  instr <- sprintf( 'rm tTestMask.nii.gz' ); system( instr )
} 
if ( file.exists( sprintf( 'tTestMask_resample.nii.gz' ) ) ) {
  instr <- sprintf( 'rm tTestMask_resample.nii.gz' ); system( instr )
} 
instr <- sprintf( '3dcalc -a brain.nii -expr \u027step(a)\u027 -prefix tTestMask.nii.gz' )
print( instr )
system( instr )
instr <- sprintf( '3dresample -inset tTestMask.nii.gz -master stats_concatenated_REML_ts_masked_Rosie.nii -prefix tTestMask_resample.nii.gz -rmode NN' )
print( instr )
system( instr )

# clean up non linear REML
if ( file.exists( sprintf( 'tOutcome_REML_nonlinear_%d.nii.gz', nBrik ) ) ) {
  instr <- sprintf( 'rm tOutcome_REML_nonlinear_%d.nii.gz', nBrik ); system( instr )
} 
instr <- paste('3dttest++',
               '-setA',
               sprintf('stats_concatenated_REML_ts_masked_Rosie.nii\u027[%d]\u027', nBrik),
               sprintf('stats_concatenated_REML_ts_masked_NWS30.nii\u027[%d]\u027', nBrik),
               sprintf('stats_concatenated_REML_ts_masked_MTR13.nii\u027[%d]\u027', nBrik),
               sprintf('stats_concatenated_REML_ts_masked_Belinda.nii\u027[%d]\u027', nBrik),
               sprintf('stats_concatenated_REML_ts_masked_Holly.nii\u027[%d]\u027', nBrik),
               sprintf('stats_concatenated_REML_ts_masked_Daniela.nii\u027[%d]\u027', nBrik),
               sprintf('stats_concatenated_REML_ts_masked_LucasEdward.nii\u027[%d]\u027', nBrik),
               sprintf('stats_concatenated_REML_ts_masked_Theo.nii\u027[%d]\u027', nBrik),
               '-mask tTestMask_resample.nii.gz',
               sprintf('-prefix tOutcome_REML_nonlinear_%d.nii.gz', nBrik)
               )
print( instr )
system( instr )

# clean up non linear 
if ( file.exists( sprintf( 'tOutcome_nonlinear_%d.nii.gz', nBrik ) ) ) {
  instr <- sprintf( 'rm tOutcome_nonlinear_%d.nii.gz', nBrik ); system( instr )
} 
instr <- paste('3dttest++',
               '-setA',
               sprintf('stats_concatenated_ts_masked_linear_Rosie.nii\u027[%d]\u027', nBrik),
               sprintf('stats_concatenated_ts_masked_linear_NWS30.nii\u027[%d]\u027', nBrik),
               sprintf('stats_concatenated_ts_masked_linear_MTR13.nii\u027[%d]\u027', nBrik),
               sprintf('stats_concatenated_ts_masked_linear_Belinda.nii\u027[%d]\u027', nBrik),
               sprintf('stats_concatenated_ts_masked_linear_Holly.nii\u027[%d]\u027', nBrik),
               sprintf('stats_concatenated_ts_masked_linear_Daniela.nii\u027[%d]\u027', nBrik),
               sprintf('stats_concatenated_REML_ts_masked_LucasEdward.nii\u027[%d]\u027', nBrik),
               sprintf('stats_concatenated_REML_ts_masked_Theo.nii\u027[%d]\u027', nBrik),
               '-mask tTestMask_resample.nii.gz',
               sprintf('-prefix tOutcome_nonlinear_%d.nii.gz', nBrik)
)
system( instr )

# clean linear REML
if ( file.exists( sprintf( 'tOutcome_REML_linear_%d.nii.gz', nBrik ) ) ) {
  instr <- sprintf( 'rm tOutcome_REML_linear_%d.nii.gz', nBrik ); system( instr )
} 
instr <- paste('3dttest++',
               '-setA',
               sprintf('stats_concatenated_REML_ts_masked_linear_Rosie.nii\u027[%d]\u027', nBrik),
               sprintf('stats_concatenated_REML_ts_masked_linear_NWS30.nii\u027[%d]\u027', nBrik),
               sprintf('stats_concatenated_REML_ts_masked_linear_MTR13.nii\u027[%d]\u027', nBrik),
               sprintf('stats_concatenated_REML_ts_masked_linear_Belinda.nii\u027[%d]\u027', nBrik),
               sprintf('stats_concatenated_REML_ts_masked_linear_Holly.nii\u027[%d]\u027', nBrik),
               sprintf('stats_concatenated_REML_ts_masked_linear_Daniela.nii\u027[%d]\u027', nBrik),
               sprintf('stats_concatenated_REML_ts_masked_LucasEdward.nii\u027[%d]\u027', nBrik),
               sprintf('stats_concatenated_REML_ts_masked_Theo.nii\u027[%d]\u027', nBrik),               
               '-mask tTestMask_resample.nii.gz',
               sprintf('-prefix tOutcome_REML_linear_%d.nii.gz', nBrik)
)
system( instr )

# clean up linear
if ( file.exists( sprintf( 'tOutcome_linear_%d.nii.gz', nBrik ) ) ) {
  instr <- sprintf( 'rm tOutcome_linear_%d.nii.gz', nBrik ); system( instr )
} 
instr <- paste('3dttest++',
               '-setA',
               sprintf('stats_concatenated_ts_masked_Rosie.nii\u027[%d]\u027', nBrik),
               sprintf('stats_concatenated_ts_masked_NWS30.nii\u027[%d]\u027', nBrik),
               sprintf('stats_concatenated_ts_masked_MTR13.nii\u027[%d]\u027', nBrik),
               sprintf('stats_concatenated_ts_masked_Belinda.nii\u027[%d]\u027', nBrik),
               sprintf('stats_concatenated_ts_masked_Holly.nii\u027[%d]\u027', nBrik),
               sprintf('stats_concatenated_ts_masked_Daniela.nii\u027[%d]\u027', nBrik),
               sprintf('stats_concatenated_REML_ts_masked_LucasEdward.nii\u027[%d]\u027', nBrik),
               sprintf('stats_concatenated_REML_ts_masked_Theo.nii\u027[%d]\u027', nBrik),
               '-mask tTestMask_resample.nii.gz',
               sprintf('-prefix tOutcome_linear_%d.nii.gz', nBrik)
)
system( instr )

}
