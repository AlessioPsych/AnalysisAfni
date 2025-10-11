rm( list=ls() ); gc();

set.seed(136);
library( plotly )

setwd('/analyse/Project0165/Fracasso_Anatomy/AHEAD_database/10007840/results/insula')

load('sub-0000_LH.RData')
storeData <- profilesArray

roiCoords <- storeData[,c(5,6,7)]
roiCoords <- t(roiCoords)
profilesArray <- storeData[,8:15]
curvatureArray <- storeData[,4]
thicknessArray <- storeData[,3]


print('........................')
print('.....check profiles.....')
print('........................')
#which( apply( abs( apply( profilesArray, 1, diff ) ), 2, sum ) < 0.00001 )
absDifferences <- abs( apply( profilesArray, 1, diff ) )
smallAbsDifferences <- apply( absDifferences < 0.0000001, 2, sum )
idxKeepProfiles <- smallAbsDifferences < 2
table( smallAbsDifferences  )

# filter profiles
print('...............................')
print('.....filter/scale profiles.....')
print('...............................')
roiCoords <- roiCoords[,idxKeepProfiles] 
curvatureArray <- curvatureArray[idxKeepProfiles] 
thicknessArray <- thicknessArray[idxKeepProfiles]
profilesArray <- profilesArray[idxKeepProfiles,]

# macro 'voxels'
print('...............................')
print('.....macro voxels..............')
print('...............................')
dataMat <- data.frame( t( roiCoords ) )
names( dataMat ) <- c('xc','yc','zc')
kMacroVoxels <- kmeans(x=dataMat, centers=500, nstart=10) 
str( dataMat )


p <- plot_ly( dataMat, x = ~xc, y = ~yc, z = ~zc, color = kMacroVoxels$cluster, size=1 ) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'xc'),
                      yaxis = list(title = 'yc'),
                      zaxis = list(title = 'zc')))
p

aggregateProfiles <- aggregate( profilesArray, by=list(kMacroVoxels$cluster), median ) # average profiles based on clustering
aggregateProfilesMat <- as.matrix( aggregateProfiles[,2:dim(aggregateProfiles)[2]] )
aggregateCoords <- aggregate( t(roiCoords), by=list(kMacroVoxels$cluster), median ) # average coordinates based on clustering

names( aggregateCoords ) <- c('cluster','xc','yc','zc')
str( aggregateCoords )

p <- plot_ly( aggregateCoords, x = ~xc, y = ~yc, z = ~zc, color = ~cluster, size=1 ) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'xc'),
                      yaxis = list(title = 'yc'),
                      zaxis = list(title = 'zc')))
p

# macro 'voxels'
print('...............................')
print('...model macro voxels..........')
print('...............................')
modelFunction <- function( idx, polyPar ) {
  singleProfile <- aggregateProfilesMat[idx, ]
  xPredictor <- seq(0,1,length.out = length(singleProfile))
  mod.lm <- lm( singleProfile ~ poly( xPredictor, polyPar ) )
  sum.mod.lm <- summary( mod.lm )
  return( as.numeric( c( coefficients( mod.lm ), sum.mod.lm$r.squared ) ) )
}
modelMacroVoxels <- sapply( 1:dim(aggregateProfilesMat)[1], modelFunction, polyPar=1 )

str( aggregateCoords )
  
p <- plot_ly( aggregateCoords, x = ~xc, y = ~yc, z = ~zc, color = modelMacroVoxels[1,], size=1 ) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'xc'),
                      yaxis = list(title = 'yc'),
                      zaxis = list(title = 'zc')))
p

keepMacroVoxel <- modelMacroVoxels[2,] < 0
modelMacroVoxels <- modelMacroVoxels[ , keepMacroVoxel ]
aggregateProfiles <- aggregateProfiles[ keepMacroVoxel, ]
aggregateProfilesMat <- aggregateProfilesMat[ keepMacroVoxel, ]
aggregateCoords <- aggregateCoords[ keepMacroVoxel, ]

print('......................................')
print('.....determine number of clusters.....')
print('......................................')
kUmap <- kmeans( aggregateProfilesMat, centers = 2, nstart = 5 )

#arrange clusters
pos01 <- apply( aggregateCoords[ kUmap$cluster==1, c(2:4) ], 2, mean )
pos02 <- apply( aggregateCoords[ kUmap$cluster==2, c(2:4) ], 2, mean )
if ( pos01[2] < pos02[2] ) {
  app <- kUmap$cluster
  kUmap$cluster[ app==1 ] <- 1
  kUmap$cluster[ app==2 ] <- 2
} else {
  app <- kUmap$cluster
  kUmap$cluster[ app==2 ] <- 1
  kUmap$cluster[ app==1 ] <- 2
}

par(mfrow=c(1,1))
p <- plot_ly( aggregateCoords, x = ~xc, y = ~yc, z = ~zc, color = kUmap$cluster, size=1, colors=c('blue','darkorange') ) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'xc'),
                      yaxis = list(title = 'yc'),
                      zaxis = list(title = 'zc')))
p

par(mfrow=c(2,2))
plot( data_umap, cex=1.5, pch=19, bty='n', xlab='umap dim. 1', ylab='umap dim. 2', las=1 )
plot( prProfiles$x[,1] ~ prProfiles$x[,2], cex=1.5, pch=19, bty='n', xlab='pr dim. 1', ylab='pr dim. 2', las=1  )
plot( data_umap, col=kUmap$cluster, cex=1.5, pch=19, bty='n', xlab='umap dim. 1', ylab='umap dim. 2', las=1 )

profile01 <- apply( aggregateProfilesMat[ kUmap$cluster==1, ], 2, mean )
profile01SE_P <- apply( aggregateProfilesMat[ kUmap$cluster==1, ], 2, mean ) + apply( aggregateProfilesMat[ kUmap$cluster==1, ], 2, sd )/2
profile01SE_M <- apply( aggregateProfilesMat[ kUmap$cluster==1, ], 2, mean ) - apply( aggregateProfilesMat[ kUmap$cluster==1, ], 2, sd )/2
profile02 <- apply( aggregateProfilesMat[ kUmap$cluster==2, ], 2, mean )
profile02SE_P <- apply( aggregateProfilesMat[ kUmap$cluster==2, ], 2, mean ) + apply( aggregateProfilesMat[ kUmap$cluster==2, ], 2, sd )/2
profile02SE_M <- apply( aggregateProfilesMat[ kUmap$cluster==2, ], 2, mean ) - apply( aggregateProfilesMat[ kUmap$cluster==2, ], 2, sd )/2
plot( profile01~seq(0,1,length.out = 8), bty='n', las=1, col='blue', xlab='cortical depth', ylab='T1w intensity (A.U.)', type='l', lwd=4,
      ylim=c( min( c( profile01SE_M, profile02SE_M ) ), max( c( profile01SE_P, profile02SE_P ) ) )  )
lines( seq(0,1,length.out = 8), profile02, col='darkorange', lwd=4 )
lines( seq(0,1,length.out = 8), profile01SE_M, col='blue', lwd=2, lty=2 )
lines( seq(0,1,length.out = 8), profile01SE_P, col='blue', lwd=2, lty=2 )
lines( seq(0,1,length.out = 8), profile02SE_M, col='darkorange', lwd=2, lty=2 )
lines( seq(0,1,length.out = 8), profile02SE_P, col='darkorange', lwd=2, lty=2 )


  
