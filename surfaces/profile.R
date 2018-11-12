args <- commandArgs(T)
print( args )

rm(list=ls())

source( sprintf('%s/AFNIio.R', '/usr/lib/afni/bin' ) )
library(quantreg)
library(boot)

#setwd('/media/alessiofracasso/storage2/prfSizeDepth_2/V4162Ben/anatomy_CBS')
setwd('/media/alessiofracasso/storage2/prfSizeDepth_2/V3814Serge/anatomy_CBS')
#setwd('/media/alessiofracasso/storage2/prfSizeDepth_2/V3813Alessio/anatomy_CBS')
#prfModel <- read.AFNI('prfModelOutput_021016_211521+orig')
prfModel <- read.AFNI('prfModelOutput_011316_174457+orig')
#prfModel <- read.AFNI('prfModelOutput_020916_195359+orig')


depth <- read.AFNI('DEPTH_al_epi_nl.nii.gz')
roi1 <- read.AFNI('LEFTV2D_al_epi_nl.nii.gz')
roi2 <- read.AFNI('RIGHTV2D_al_epi_nl.nii.gz')
#depth <- read.AFNI('DEPTH_al_epi.nii.gz')
#roi1 <- read.AFNI('LEFTV1_al_epi.nii.gz')
#roi2 <- read.AFNI('RIGHTV1_al_epi.nii.gz')

selVoxels <- function(prfModel, index, roi1, roi2) {
  vol <- prfModel$brk[,,,index]
  out <- vol[ roi1$brk | roi2$brk == 1 ]
  return(out)
}

area <- sapply( seq(1,8), selVoxels, prfModel=prfModel, roi1=roi1, roi2=roi2  )
selDepth <- selVoxels(depth,1,roi1,roi2)
areaFull <- data.frame( cbind(area, selDepth) )
names( areaFull ) <- c('coherence', 'sigmaUpCenter', 'sigmaUpSurround', 
                       'sigma', 'ecc', 'betaIntercept', 'betaSlope1', 
                       'betaSlope2', 'distance')
  
flagData <- areaFull$coherence > quantile(areaFull$coherence, 0.3 ) & areaFull$coherence <= quantile(areaFull$coherence, 0.9 ) &
    areaFull$sigmaUpCenter<15 & areaFull$sigmaUpSurround<15 &
    areaFull$sigmaUpCenter>0.2 & areaFull$sigmaUpSurround>0.2 &
    areaFull$ecc > 0.5 & areaFull$ecc < 5 &
    areaFull$betaSlope1 > 0 & areaFull$betaSlope1 < 30 &
    areaFull$betaSlope2 > 0 & areaFull$betaSlope2 < 30 &
    areaFull$distance>0.05 & areaFull$distance<0.95

areaFilt <- areaFull[ flagData, ]

x11()
nData <- names(areaFilt)
par(mfrow=c(3,4))
for (k in 1:9) {
  hist( areaFilt[,k], main=nData[k] )
  if (k>9) { par(new=TRUE) }
}
plot( areaFilt$sigmaUpCenter ~ areaFilt$ecc, cex=0.1 )
abline( rq( areaFilt$sigmaUpCenter ~ areaFilt$ecc, tau=0.5 ), lwd=2, col='red'  )
plot( areaFilt$sigmaUpSurround ~ areaFilt$ecc, cex=0.1 )
abline( rq( areaFilt$sigmaUpSurround ~ areaFilt$ecc, tau=0.5 ), lwd=2, col='red'  )

index <- seq(1:dim(areaFilt)[1])
depthCoeff <- function( areaFilt, index, pDistance ) {    
  dataFilt <- areaFilt[index, ]
  distanceCut <- as.numeric( cut( dataFilt$distance, quantile( dataFilt$distance, probs=seq(0,1,pDistance) ), include.lowest=TRUE ) )  
  lmDepth <- function( dataFilt, distanceCut, indexDepth, indexSigma ) {
    x <- dataFilt$ecc[distanceCut==indexDepth]
    y <- dataFilt[ distanceCut==indexDepth, indexSigma ]    
    xCut <- as.numeric( cut( x, unique( quantile( x, probs=seq(0,1,0.1) ) ), include.lowest=TRUE ) )
    xEst <- array( tapply( x, list(xCut), mean ) )
    yEst <- array( tapply( y, list(xCut), mean ) )
    mod <- summary( lm( yEst ~ xEst ) )
    coeff <- coefficients(mod)[,1]
    #mod <- summary( ( y ~ x ) )
    #coeff <- coefficients(mod)[,1]
    return( coeff )
  } 
  counter <- seq( 1, length( seq(0,1,pDistance) )-1 )
  centerPar <- sapply( counter, lmDepth, dataFilt=dataFilt, distanceCut=distanceCut, indexSigma=2 )
  surrPar <- sapply( counter, lmDepth, dataFilt=dataFilt, distanceCut=distanceCut, indexSigma=3 )
  #return( ( rbind( centerPar, surrPar ) ) )
  return( round( ( c( centerPar[1,], centerPar[2,], surrPar[1,], surrPar[2,] ) ), 3 ) )
}

pDistance <- 0.1
nElementsProfile <- length( seq(0,1,pDistance) ) - 1
a <- depthCoeff( areaFilt, index, pDistance)


distanceCut <- as.numeric( cut( areaFilt$distance, quantile( areaFilt$distance, probs=seq(0,1,pDistance) ), include.lowest=TRUE ) )  
nRep <- 100
bootCoeff <- array(0,c(4,nElementsProfile,nRep))
bootMod <- boot( areaFilt, depthCoeff, pDistance=pDistance, nRep, strata=distanceCut)
for (k in 1:nRep ) {
  bootCoeff[,,k] <- t( matrix(bootMod$t[k,], ncol=4, nrow=nElementsProfile ) )
}
multFactor <- seq(2.9,3.1,0.1)
center <- array(0,dim(bootCoeff)[2:3])
surround <- array(0,dim(bootCoeff)[2:3])
for (n in 1:length(multFactor)) {
  center <- center + bootCoeff[1,,] + bootCoeff[2,,]*multFactor[n]
  surround <- surround + bootCoeff[3,,] + bootCoeff[4,,]*multFactor[n]
}
center <- center / length(multFactor)
surround <- surround / length(multFactor)

xVar <- tapply( areaFilt$distance, list(distanceCut), median )
plot( apply( center, 1, quantile, probs=0.5 ) ~ xVar  )
lines( xVar, apply( center, 1, quantile, probs=0.975 ), lwd=1, lty=2  )
lines( xVar, apply( center, 1, quantile, probs=0.025 ), lwd=1, lty=2  )
spCenter <- smooth.spline( xVar, apply( center, 1, quantile, probs=0.5 ), spar=0.4 )
lines(spCenter)

plot( apply( surround, 1, quantile, probs=0.5 ) ~ xVar )
lines( xVar, apply( surround, 1, quantile, probs=0.975 ), lwd=1, lty=2  )
lines( xVar, apply( surround, 1, quantile, probs=0.025 ), lwd=1, lty=2  )
spSurround <- smooth.spline( xVar, apply( surround, 1, quantile, probs=0.5 ), spar=0.4 )
lines(spSurround)

lTrend <- seq(-1,1,length.out=10)
qTrend <- lTrend^2

summary( lm( spCenter$y ~ lTrend + qTrend ) )
summary( lm( spSurround$y ~ lTrend + qTrend ) )


distanceCut <- cut( areaFilt$distance, quantile( areaFilt$distance, probs=seq(0,1,0.1) ) )
plot( tapply( areaFilt$betaSlope1, list(distanceCut), mean ) )
plot( tapply( areaFilt$coherence, list(distanceCut), mean ) )
plot( tapply( areaFilt$sigmaUpCenter, list(distanceCut), median ) )
plot( tapply( areaFilt$sigmaUpSurround, list(distanceCut), median ) )