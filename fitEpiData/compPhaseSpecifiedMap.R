args <- commandArgs(T)
print( args )

#rm(list=ls())
#setwd('/home/fracasso/data/stripesPSF/participants/S01')
#args <- c('trWave.nii.gz', '3.14', '/packages/afni/16.1.13','/home/fracasso/analysisAfni/surfaces')

source( sprintf('%s/AFNIio.R', args[3] ) )
source( sprintf('%s/coordinateFromLinearIndex.r', args[4] ) )
source( sprintf('%s/linearIndexFromCoordinate.r', args[4] ) )
instr <- sprintf('3dcopy %s ttt_data.nii.gz', args[1])
system( instr )
refPhase <- args[2]

instr <- sprintf( '3dcalc -a ttt_data.nii.gz[0] -b ttt_data.nii.gz[1] -c ttt_data.nii.gz[2] -expr \u0027a*cos(c-%s)\u0027 -prefix ttt_amp.nii.gz', refPhase )
system( instr )
instr <- sprintf( '3dcalc -a ttt_data.nii.gz[0] -b ttt_data.nii.gz[1] -c ttt_data.nii.gz[2] -expr \u0027b*cos(c-%s)\u0027 -prefix ttt_co.nii.gz', refPhase )
system( instr )
instr <- '3dTcat -prefix phaseSpecCoh.nii.gz ttt_co.nii.gz ttt_amp.nii.gz'
system( instr )
instr <- '3drefit -sublabel 0 coh phaseSpecCoh.nii.gz'
system( instr )
instr <- '3drefit -sublabel 1 amp phaseSpecCoh.nii.gz'
system( instr )

system('rm ttt_*')

