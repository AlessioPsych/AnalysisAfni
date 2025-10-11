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
#profilesArray <- profilesArray
profilesArray <- sapply( 1:10, FUN = cleanData, simplify = TRUE  )

aggregateCoords_01 <- dfOutSel[,c(3,4,5)]

# Determine number of clusters and some plots
print('......................................')
print('.....determine number of clusters.....')
print('......................................')

#kUmap <- kmeans( aggregateProfilesMat_01, centers = 2, nstart = 25,  algorithm = c("Hartigan-Wong" ), iter.max=1500 )
kUmap <- kmeans( profilesArray, centers = 2, nstart = 25,  algorithm = c("Hartigan-Wong" ), iter.max=1500 )

#arrange clusters
pos01 <- apply( aggregateCoords_01[ kUmap$cluster==1, c(1:3) ], 2, mean )
pos02 <- apply( aggregateCoords_01[ kUmap$cluster==2, c(1:3) ], 2, mean )
if ( pos01[2] < pos02[2] ) {
  app <- kUmap$cluster
  kUmap$cluster[ app==1 ] <- 1
  kUmap$cluster[ app==2 ] <- 2
} else {
  app <- kUmap$cluster
  kUmap$cluster[ app==2 ] <- 1
  kUmap$cluster[ app==1 ] <- 2
}

p <- plot_ly( aggregateCoords_01, x = ~i, y = ~j, z = ~k, color = kUmap$cluster, size=1, colors=c('blue','darkorange') ) %>%
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