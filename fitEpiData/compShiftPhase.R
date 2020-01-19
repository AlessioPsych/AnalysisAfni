args <- commandArgs(T)
print( args )

#rm(list=ls())
#setwd('/analyse/Project0226/dataToSort/eyeMovPoster/FSH19_KastnerVest/kastner_ts_mc')
#args <- c('trWave.nii.gz', '30', 'trWave_shift01.nii.gz' )

AFNI_INSTALLDIR <- Sys.getenv( x = 'AFNI_INSTALLDIR' )
generalPurpose_DIR <- Sys.getenv( x = 'AFNI_TOOLBOXDIRGENERALPURPOSE' )
surfaces_Dir <- Sys.getenv( x = 'AFNI_TOOLBOXDIRSURFACES' )
source( sprintf('%s/AFNIio.R', AFNI_INSTALLDIR ) )
source( sprintf('%s/coordinateFromLinearIndex.r', surfaces_Dir ) )
source( sprintf('%s/linearIndexFromCoordinate.r', surfaces_Dir ) )
source( sprintf('%s/scaleData.R', generalPurpose_DIR ) )
library( pracma )

# arrange inputs
trInput <- args[1]
shiftInput <- as.numeric( args[2] )
trOutput <- args[3]

# load inputs
shiftRadiants <- deg2rad( shiftInput )
dataFile <- read.AFNI( filename=trInput )
dataAmp <- dataFile$brk[,,,1]
dataCo <- dataFile$brk[,,,2]
dataPh <- dataFile$brk[,,,3]

# test modulus:
#( c( 6.28, 6.2, 3 ) + 0.1 ) %% (2*pi)
#( c( 0, 3.24, 0.5 ) - 0.1 ) %% (2*pi)

# shift phase backwards: positive shift value, to account for hemodynamic delay
# phi <- deg2rad( +30 )
# piSpace <- seq(0,2*pi,0.2)
# plot( sin( (piSpace+phi) ) ~ piSpace, type='l', col='red'  ); lines( piSpace, sin( piSpace ) ); legend(4,0.5, legend=c('original','shifted'), fill=c('black','red')  )

# shift phase forwards: negative shift value, does not account for hemodynamic delay
# phi <- deg2rad( -30 )
# piSpace <- seq(0,2*pi,0.2)
# plot( sin( (piSpace+phi) ) ~ piSpace, type='l', col='red'  ); lines( piSpace, sin( piSpace ) ); legend(4,0.5, legend=c('original','shifted'), fill=c('black','red')  )

# shift phase data
dataPhShift <- ( dataPh + shiftRadiants ) %% (2*pi)
dataFile$brk[,,,3] <- dataPhShift

x11( width=6, height=3 )
par(mfrow=c(1,3))
thrValue <- 0.25 #for visualization only
hist( dataPh[ dataCo>thrValue ], xlim=c(0,9), xlab='original ph', main=sprintf('co thr = %1.2f', thrValue) )
hist( ( dataPh[ dataCo>thrValue ] + shiftRadiants ), xlim=c(0,9), xlab='original ph, shifted, no modulus', main=sprintf('co thr = %1.2f', thrValue) )
hist( dataPhShift[ dataCo>thrValue ], xlim=c(0,9), xlab='original ph, shifted', main=sprintf('co thr = %1.2f', thrValue) )

print('prepare and save output...')

volFileName <- sprintf( trOutput )
write.AFNI(volFileName,
           brk=dataFile$brk,
           view='+orig',
           orient=dataFile$orient,
           origin=dataFile$origin,
           defhead=dataFile$NI_head )


