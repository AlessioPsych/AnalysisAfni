args <- commandArgs(T)
print( args )

#setwd('/home/fracasso/data/pRF/sub1')
#args <- c('prfParameters.1D')

prfFile <- args[1]

library( abind )
library( pracma )

fileData <- scan( prfFile, what=list('') )

baselineLength <- as.numeric( fileData[[1]][2] )
baselineLengthBegin <- as.numeric( fileData[[1]][4] )
baselineLengthEnd <- as.numeric( fileData[[1]][6] )
bWidth <- as.numeric( fileData[[1]][8] )
bStep <- as.numeric( fileData[[1]][10] )
nSteps <- as.numeric( fileData[[1]][12] )
width <- as.numeric( fileData[[1]][14] )

# build stimuli sequence
print('build stimuli... (1)')

x <- seq(-width/2,width/2,length.out = 100 )
y <- x
bStart <- -width/2
stimHorizontal <- array( 0, c( length(x), length(y), nSteps ) )
mesh2dCoords <- meshgrid( x, y ) # build the mask
circleMask <- mesh2dCoords$X^2 + mesh2dCoords$Y^2 < (width/2)^2 # mask
for (k in 1:nSteps) {
  if (k == 1) {
    xStart <- bStart
    xEnd <- xStart + bWidth
    xSel <- x>=xStart & x<xEnd
  }
  if (k>1) {
    xStart <- xStart + bStep
    xEnd <- xEnd + bStep
    xSel <- x>=xStart & x<xEnd    
  }
  stimHorizontal[ , x>=xStart & x<xEnd, k] <- 1 
  stimHorizontal[,,k] <- stimHorizontal[,,k]*circleMask
}

print('build stimuli... (2)')
phi <- deg2rad(-45)
mesh2dCoords <- meshgrid( x, y )
xArray <- array( mesh2dCoords$X )
yArray <- array( mesh2dCoords$Y )
matTransform <- matrix( c( cos(phi), sin(phi), -sin(phi), cos(phi) ), byrow=TRUE, ncol=2 )
coordsTransformedMin45 <- matTransform%*%rbind(xArray,yArray)

phi <- deg2rad(45)
mesh2dCoords <- meshgrid( x, y )
xArray <- array( mesh2dCoords$X )
yArray <- array( mesh2dCoords$Y )
matTransform <- matrix( c( cos(phi), sin(phi), -sin(phi), cos(phi) ), byrow=TRUE, ncol=2 )
coordsTransformedPlus45 <- matTransform%*%rbind(xArray,yArray)

stimAngleMin45 <- array( 0, c( length(x), length(y), nSteps ) )
stimAnglePlus45 <- array( 0, c( length(x), length(y), nSteps ) )
for (k in 1:nSteps) {
  img <- stimHorizontal[,,k]
  
  img02 <- interp2( x, y, img, coordsTransformedMin45[1,], coordsTransformedMin45[2,], method=c('nearest')  )
  img03 <- array( img02, dim( stimHorizontal[,,1] ) )
  img03[ is.na(img03) ] <- 0
  stimAngleMin45[,,k] <- img03
  
  img02 <- interp2( x, y, img, coordsTransformedPlus45[1,], coordsTransformedPlus45[2,], method=c('nearest')  )
  img03 <- array( img02, dim( stimHorizontal[,,1] ) )
  img03[ is.na(img03) ] <- 0
  stimAnglePlus45[,,k] <- img03
  
  img02 <- interp2( x, y, img, xArray, yArray, method=c('nearest')  )
  img03 <- array( img02, dim( stimHorizontal[,,1] ) )
  img03[ is.na(img03) ] <- 0
  stimHorizontal[,,k] <- img03
}

print('build stimuli... (3)')
stimUp <- stimHorizontal
stimDown <- stimHorizontal[,,nSteps:1]
stimLeft <- aperm( stimHorizontal, c(2,1,3) )
stimRight <- stimLeft[,,nSteps:1]
stimUpAngleDx <- stimAngleMin45
stimUpAngleSx <- stimAnglePlus45
stimDownAngleDx <- stimAnglePlus45[,,nSteps:1]
stimDownAngleSx <- stimAngleMin45[,,nSteps:1]

print('build stimuli... (4)')
baselineStep <- array( 0, c( length(x), length(y), baselineLength ) )
baselineStepBegin <- array( 0, c( length(x), length(y), baselineLengthBegin ) )
baselineStepEnd <- array( 0, c( length(x), length(y), baselineLengthEnd ) )

instrBeg <- 'stimSeq <- abind( '
indexElements <- 15:length( fileData[[1]] )
for (k in indexElements) {
  if (k==15) {
    instrPart <- sprintf( ' %s ,', fileData[[1]][ k ] )
  }
  if (k>15 & k<max(indexElements)) {
    instrPart <- strcat( instrPart, sprintf( '%s ,', fileData[[1]][ k ] ) )
  }
  if (k==max(indexElements)) {
    instrPart <- strcat( instrPart, sprintf( '%s )', fileData[[1]][ k ] ) )
  }
}
instr <- strcat( instrBeg, instrPart ) 
eval(parse(text=instr))

print('save stimuli... ')
save( stimSeq, x, y, file='prfStimuli.RData' )
