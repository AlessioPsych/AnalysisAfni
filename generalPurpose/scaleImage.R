args <- commandArgs(T)
print( args )

#debug: 
# setwd('/media/alessiofracasso/storage2/dataLaminar_leftVSrightNEW/V2677leftright_Ric/resultsCBS/coregistration')
# args <- c('coreg_tempFile_MPRAGE_shft_al+orig', '/home/alessiofracasso/Dropbox/analysisAfni', '/usr/lib/afni/bin', 'invertedContrast.nii.gz' )

eVar <- Sys.getenv(c('AFNI_TOOLBOXDIR','AFNI_INSTALLDIR'))

source( sprintf('%s/scaleData.R', args[5] ) )
source( sprintf('%s/AFNIio.R', args[6] ) )

volume <- read.AFNI( args[1] )

volumeData <- volume$brk[,,,1]
maskData <- volumeData!=0

oldMin <- min( volumeData[maskData] )
oldMax <- max( volumeData[maskData] )

newMax <- as.numeric( args[2] )
newMin <- as.numeric( args[3] )


dataInverted <- scaleData( volumeData[maskData], newMax, newMin )
volumeDataInverted <- array( 0, dim( volumeData ) )
volumeDataInverted[maskData] <- dataInverted

write.AFNI(args[4],
           brk=volumeDataInverted,
           label=NULL,
           view='+orig',
           orient=volume$orient,
           origin=volume$origin,
           delta=volume$delta,
           defhead=volume$NI_head )
