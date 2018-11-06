args <- commandArgs(T)
print( args )
library(pracma)

#debug: 
#setwd('/data1/projects/UMC_data_Fracasso/tests/MP2RAGE_test_Gla/wetransfer-57c266/sub-SLW06_T1_ISO_04mm')
#args <- c('_levelset_WM_levelset.nii.gz','probFile.nii.gz','4','/packages/afni/17.0.13')

#source( sprintf('%s/scaleData.R', args[5] ) )
source( sprintf('%s/AFNIio.R', args[4] ) )

volume <- read.AFNI( args[1] )
volumeBrk_prob <- sigmoid( volume$brk[,,,1]*-1, as.numeric( args[3] ), 0 )

write.AFNI(args[2],
           brk=volumeBrk_prob,
           label=NULL,
           view='+orig',
           orient=volume$orient,
           origin=volume$origin,
           delta=volume$delta,
           defhead=volume$NI_head )

# set zero to 0.00001


#plot( sigmoid( seq(-5,5,0.1), 1, 0 )~seq(-5,5,0.1) )