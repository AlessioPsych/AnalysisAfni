args <- commandArgs(T)
print( args )

#rm(list=ls())
#setwd('/home/fracasso/Desktop/coregistration')
#args <- c('warpTop_PLUS+orig','t1_mask_no_blur.nii.gz','/packages/afni/16.1.13')

source( sprintf('%s/AFNIio.R', args[3] ) )
targetFilename <- args[1]
masterFilename <- args[2]

instr <- sprintf('3dinfo -prefix_noext %s > tttemporary.1D', targetFilename )
system( instr )
cleanName <- scan( 'tttemporary.1D', n=1, what=list('character') )[[1]]
instr <- sprintf( '3dcopy %s %s_refit.nii.gz', args[1], cleanName )
system(instr)

masterFile <- read.AFNI( masterFilename  )
chVect <- as.character(masterFile$NI_head$IJK_TO_DICOM_REAL$dat)
instr <- sprintf('3drefit -atrfloat IJK_TO_DICOM_REAL \u027%s %s %s %s %s %s %s %s %s %s %s %s\u0027 %s_refit.nii.gz', chVect[1],chVect[2],chVect[3],chVect[4],chVect[5],chVect[6],chVect[7],chVect[8],chVect[9],chVect[10],chVect[11],chVect[12],cleanName)
system( instr )

system('rm tttemporary.1D')