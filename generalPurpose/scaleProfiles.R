args <- commandArgs(T)
print( args )

#debug: 
# setwd('/media/alessiof/Data/anatTestCode')
# args <- c('profilesOut_profiles.nii.gz', '8000', '1000', 'profilesScaled01.nii.gz' )

eVar <- Sys.getenv(c('AFNI_TOOLBOXDIR','AFNI_INSTALLDIR','AFNI_TOOLBOXDIRGENERALPURPOSE'))

source( sprintf('%s/scaleData.R', eVar[3] ) )
source( sprintf('%s/AFNIio.R', eVar[2] ) )

print('...........................')
print('...........................')
print('...........................')
print('reading input volume.......')
print('...........................')
print('...........................')
print('...........................')
volume <- read.AFNI( args[1] )

print('............................')
print('............................')
print('keep valid voxels and.......')
print('convert.....................')
print('input volume to array.......')
print('............................')
print('............................')
volumeData <- volume$brk[,,,1]
maskData <- volumeData>0.00001
tempProfiles <- array( 999, c( sum(maskData), dim(volume$brk)[4] ) )
for (nVols in 1:dim(volume$brk)[4] ) {
  tempVol <- volume$brk[,,,nVols]
  tempProfiles[,nVols] <- tempVol[ maskData ]
}

print('............................')
print('............................')
print('scale valid voxels..........')
print('............................')
print('............................')
newMax <- as.numeric( args[2] )
newMin <- as.numeric( args[3] )
dataScaled <- scaleData( tempProfiles, newMax, newMin )

print('............................')
print('............................')
print('arrange output volume.......')
print('save........................')
print('............................')
volumeDataOut <- array( 0, dim( volume$brk ) )
for (nVols in 1:dim(volume$brk)[4] ) {
  tempVolArray <- dataScaled[,nVols]
  tempVolVolume <- array( 0, dim( volumeDataOut )[1:3]  )
  tempVolVolume[ maskData ] <- tempVolArray
  volumeDataOut[,,,nVols ] <- tempVolVolume
}

write.AFNI(args[4],
           brk=volumeDataOut,
           label=NULL,
           view='+orig',
           orient=volume$orient,
           origin=volume$origin,
           delta=volume$delta,
           defhead=volume$NI_head )
