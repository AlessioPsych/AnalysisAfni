rm( list=ls() ); gc();
graphics.off();


instr <- 'Rscript 09_insula_test_analysis_hemiCorrection_vonEconomo_Part1.R'
system( instr )
rm( list=ls() ); gc();
graphics.off();
instr <- 'Rscript 09_insula_test_analysis_hemiCorrection_vonEconomo_Part2.R'
system( instr )
