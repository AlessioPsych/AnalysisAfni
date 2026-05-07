rm( list=ls() )
set.seed(1);
afniDir <- Sys.getenv(x='AFNI_INSTALLDIR')
afniSurfaces <- Sys.getenv(x='AFNI_TOOLBOXDIRSURFACES')
source( sprintf( '%s/AFNIio.R', afniDir ) )
source( sprintf( '%s/coordinateFromLinearIndex.r', afniSurfaces ) )
library( plotly )
library( ica )
library( dbscan )

#saveRoi <- '/analyse/Project0165/Fracasso_Anatomy/AHEAD_database/10007840/results/insula'
#appendFilename <- 'LH'

mainDir <- '/analyse/Project0165/Fracasso_Anatomy/AHEAD_database/10007840'
partDir <- 'part1/'

setwd( mainDir )
setwd( partDir )
getwd()
participantsFolders <- dir()

#for ( participantsFoldersIndex in 2:4 ) { #seq( 1,length( participantsFolders )

participantsFoldersIndex <- 2
setwd( mainDir )
setwd( partDir )
print( '...' )
print( '...' )
print( sprintf( '%s...', participantsFolders[ participantsFoldersIndex ] ) )
print( '...' )
print( '...' )
setwd( sprintf( '%s', participantsFolders[ participantsFoldersIndex ] ) )
setwd( dir()[1] )
setwd( dir()[1] )
print( dir() )
setwd( 'AHEAD_test' )
setwd( 'SUMA' )
print( dir() )

atlas_lh <- read.table('std.141.aparc+aseg_REN_all_midpoint.lh.1D.dset', header=FALSE)
names( atlas_lh ) <- c('node', 'index', 'i', 'j', 'k', 'vals', 'atlasIdx' )

curvature_lh <- read.table('std.141.lh.curv_ppp_.1D.dset', header=FALSE)
names( curvature_lh ) <- c('node', 'curvature' )

sulc_lh <- read.table('std.141.lh.sulc_ppp_.1D.dset', header=FALSE)
names( sulc_lh ) <- c('node', 'sulc' )

thickness_lh <- read.table('std.141.lh.thickness_ppp_.1D.dset', header=FALSE)
names( thickness_lh ) <- c('node', 'thickness' )

t1Profiles_lh <- read.table( 'std.141.t1w_segvals_10.lh.1D.dset', header=FALSE )
names( t1Profiles_lh ) <- c('node', 'index', 'i', 'j', 'k', 'vals', 'V0', 'V1', 'V2', 'V3', 'V4', 'V5', 'V6', 'V7', 'V8', 'V9' )
str( t1Profiles_lh )

insulaLeft <- data.frame( t1Profiles_lh[ atlas_lh$atlasIdx == 80, ],
                          curvature_lh$curvature[ atlas_lh$atlasIdx == 80 ], 
                          thickness_lh$thickness[ atlas_lh$atlasIdx == 80 ],
                          sulc_lh$sulc[ atlas_lh$atlasIdx == 80 ]) 
names( insulaLeft ) <- c('node', 'index', 'i', 'j', 'k', 'vals', 'V0', 'V1', 'V2', 'V3', 'V4', 'V5', 'V6', 'V7', 'V8', 'V9', 'curvature', 'thickness', 'sulc' )
str( insulaLeft )

rm(list=ls())
setwd('/analyse/Project0165/Fracasso_Anatomy/AHEAD_database/10007840/results/insulaSurface')
load( 'sub-0083.RData' )

dfOutSel <- dataOut[ dataOut$modality=='t1w' & dataOut$hemi=='LH',  ]

profilesArray <- dfOutSel[,7:16]
roiCoords <- t( array( as.matrix( dfOutSel[,3:5] ), c( dim(dfOutSel)[1], 3 ) ) )

cleanData <- function( idxIn ) {
  #idxIn <- 4
  modCovar <- lm( profilesArray[,idxIn] ~ poly( dfOutSel$thickness, 2 ) * poly( dfOutSel$curvature, 2 ) )
  summary( modCovar )
  dataClean <- round( modCovar$resid + coefficients( modCovar )[1], 2 )
  #dataClean <- round( modCovar$resid, 2 )
  return( dataClean )
}
profilesArray <- profilesArray
profilesArray <- sapply( 1:10, FUN = cleanData, simplify = TRUE  )

# check profiles:
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
profilesArray <- profilesArray[idxKeepProfiles,]
dfOutSel <- dfOutSel[idxKeepProfiles,] 

# macro 'voxels'
print('...............................')
print('.....macro voxels..............')
print('...............................')
dataMat <- data.frame( t( roiCoords ) )
names( dataMat ) <- c('xc','yc','zc')
kMacroVoxels <- kmeans(x=dataMat, centers=2000, nstart=10) 
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
aggregateCurvature <- aggregate( dfOutSel$curvature, by=list(kMacroVoxels$cluster), median ) # average curvature based on clustering
aggregateThickness <- aggregate( dfOutSel$thickness, by=list(kMacroVoxels$cluster), median ) # average curvature based on clustering
aggregateCoords$curvature <- aggregateCurvature[,2]
aggregateCoords$thickness <- aggregateThickness[,2]

names( aggregateCoords ) <- c('cluster','xc','yc','zc', 'curvature', 'thickness')
str( aggregateCoords )
p <- plot_ly( aggregateCoords, x = ~xc, y = ~yc, z = ~zc, color = ~cluster, size=1 ) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'xc'),
                      yaxis = list(title = 'yc'),
                      zaxis = list(title = 'zc')))
p

p <- plot_ly( aggregateCoords, x = ~xc, y = ~yc, z = ~zc, color = ~curvature, size=1 ) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'xc'),
                      yaxis = list(title = 'yc'),
                      zaxis = list(title = 'zc')))
p

p <- plot_ly( aggregateCoords, x = ~xc, y = ~yc, z = ~zc, color = ~thickness, size=1 ) %>%
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
modelMacroVoxels <- sapply( 1:dim(aggregateProfilesMat)[1], modelFunction, polyPar=2 )

str( aggregateCoords )
p <- plot_ly( aggregateCoords, x = ~xc, y = ~yc, z = ~zc, color = modelMacroVoxels[3,], size=1 ) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'xc'),
                      yaxis = list(title = 'yc'),
                      zaxis = list(title = 'zc')))
p

colorMacroVoxel <- as.numeric( modelMacroVoxels[2,] < 0 )
p <- plot_ly( aggregateCoords, x = ~xc, y = ~yc, z = ~zc, color = colorMacroVoxel, size=1 ) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'xc'),
                      yaxis = list(title = 'yc'),
                      zaxis = list(title = 'zc')))
p

#####################################################################
keepMacroVoxel <-  modelMacroVoxels[2,] < 0 # rep( TRUE, length( modelMacroVoxels[2,] ) ) #
modelMacroVoxels_01 <- modelMacroVoxels[ , keepMacroVoxel ]
aggregateProfiles_01 <- aggregateProfiles[ keepMacroVoxel, ]
aggregateProfilesMat_01 <- aggregateProfilesMat[ keepMacroVoxel, ]
aggregateCoords_01 <- aggregateCoords[ keepMacroVoxel,  ]
#####################################################################

#kNNdistplot( aggregateProfilesMat_01, k=1 )
#abline( h=270, lwd=2, lty=2 )
#dbscan( aggregateProfilesMat_01, 100, 10 )

pc <- prcomp( aggregateProfilesMat_01, center=FALSE, scale=FALSE )
summary( pc )

p <- plot_ly( aggregateCoords_01, x = ~xc, y = ~yc, z = ~zc, size=1 ) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'xc'),
                      yaxis = list(title = 'yc'),
                      zaxis = list(title = 'zc')))
p

x11()
par(mfrow=c(2,2))
plot( modelMacroVoxels_01[2,] ~ modelMacroVoxels_01[1,], bty='n', las=1, xlab='intercept', ylab='slope', cex=1, pch=19, col=rgb(0.5,0.5,0.5,0.5) )
abline( lm( modelMacroVoxels_01[2,] ~ modelMacroVoxels_01[1,] ), lwd=2, col='red' )
hist( modelMacroVoxels_01[1,] )  
hist( modelMacroVoxels_01[2,] )  

# Determine number of clusters and some plots
print('......................................')
print('.....determine number of clusters.....')
print('......................................')

#kUmap <- kmeans( aggregateProfilesMat_01, centers = 2, nstart = 25,  algorithm = c("Hartigan-Wong" ), iter.max=1500 )
kUmap <- kmeans( aggregateProfilesMat_01, centers = 2, nstart = 25,  algorithm = c("Hartigan-Wong" ), iter.max=1500 )

#arrange clusters
pos01 <- apply( aggregateCoords_01[ kUmap$cluster==1, c(2:4) ], 2, mean )
pos02 <- apply( aggregateCoords_01[ kUmap$cluster==2, c(2:4) ], 2, mean )
if ( pos01[2] < pos02[2] ) {
  app <- kUmap$cluster
  kUmap$cluster[ app==1 ] <- 1
  kUmap$cluster[ app==2 ] <- 2
} else {
  app <- kUmap$cluster
  kUmap$cluster[ app==2 ] <- 1
  kUmap$cluster[ app==1 ] <- 2
}

p <- plot_ly( aggregateCoords_01, x = ~xc, y = ~yc, z = ~zc, color = kUmap$cluster, size=1, colors=c('blue','darkorange') ) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'xc'),
                      yaxis = list(title = 'yc'),
                      zaxis = list(title = 'zc')))
p

x11( width=4, height=4 )
par(mfrow=c(1,1))
profile01 <- apply( aggregateProfilesMat_01[ kUmap$cluster==1, ], 2, mean )
profile01SE_P <- apply( aggregateProfilesMat_01[ kUmap$cluster==1, ], 2, mean ) + apply( aggregateProfilesMat_01[ kUmap$cluster==1, ], 2, sd )/2
profile01SE_M <- apply( aggregateProfilesMat_01[ kUmap$cluster==1, ], 2, mean ) - apply( aggregateProfilesMat_01[ kUmap$cluster==1, ], 2, sd )/2
profile02 <- apply( aggregateProfilesMat_01[ kUmap$cluster==2, ], 2, mean )
profile02SE_P <- apply( aggregateProfilesMat_01[ kUmap$cluster==2, ], 2, mean ) + apply( aggregateProfilesMat_01[ kUmap$cluster==2, ], 2, sd )/2
profile02SE_M <- apply( aggregateProfilesMat_01[ kUmap$cluster==2, ], 2, mean ) - apply( aggregateProfilesMat_01[ kUmap$cluster==2, ], 2, sd )/2
plot( profile01~seq(0,1,length.out = 10), bty='n', las=1, col='blue', xlab='cortical depth', ylab='T1w intensity (A.U.)', type='l', lwd=4,
      ylim=c( min( c( profile01SE_M, profile02SE_M ) ), max( c( profile01SE_P, profile02SE_P ) ) ), axes=FALSE, cex.lab=1.5  )
axis(1,seq(0,1,0.5), cex.axis=1.5); axis(2, seq( round(min( c( profile01SE_M, profile02SE_M ) ),0), round(max( c( profile01SE_P, profile02SE_P ) ),0), length.out = 3 ), cex.axis=1.5 )
lines( seq(0,1,length.out = 10), profile02, col='darkorange', lwd=4 )
lines( seq(0,1,length.out = 10), profile01SE_M, col='blue', lwd=2, lty=2 )
lines( seq(0,1,length.out = 10), profile01SE_P, col='blue', lwd=2, lty=2 )
lines( seq(0,1,length.out = 10), profile02SE_M, col='darkorange', lwd=2, lty=2 )
lines( seq(0,1,length.out = 10), profile02SE_P, col='darkorange', lwd=2, lty=2 )





#}