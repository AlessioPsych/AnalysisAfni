rm( list=ls() )

instr <- 'Rscript 05_extractRoiData_insula_LH_part01_vonEconomo.R'
system( instr )
instr <- 'Rscript 06_extractRoiData_insula_RH_part01_vonEconomo.R'
system( instr )
instr <- 'Rscript 07_extractRoiData_insula_LH_part02_vonEconomo.R'
system( instr )
instr <- 'Rscript 08_extractRoiData_insula_RH_part02_vonEconomo.R'
system( instr )