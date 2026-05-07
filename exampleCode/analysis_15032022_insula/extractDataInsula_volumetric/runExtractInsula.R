rm( list=ls() )

instr <- 'Rscript 05_extractRoiData_insula_LH_part01.R'
system( instr )
instr <- 'Rscript 06_extractRoiData_insula_RH_part01.R'
system( instr )
instr <- 'Rscript 07_extractRoiData_insula_LH_part02.R'
system( instr )
instr <- 'Rscript 08_extractRoiData_insula_RH_part02.R'
system( instr )