args <- commandArgs(T)
print( args )

#setwd('/analyse/Project0226/dataRepository/anatomies/surfaceAtlases/suma_MNI152_2009_princetonAtlas')
#args <- c('0.5', 'testRot','7_0.4_axis_lh_axis.dset')

# input
rot <- as.numeric( args[1] )
outName <- args[2]
inputFile <- args[3]

# get data
vertexF_in <- read.table( inputFile, header=FALSE )
nodes <- vertexF_in[,1]
vertexF_in <- t( vertexF_in[,c(2,3)] )
n_in <- dim( vertexF_in )[2]

# center and rotate
idxCenter <- which.min( abs( vertexF_in[1,] )  )
vertexFC <- vertexF_in - t( array( rep( vertexF_in[ , idxCenter ], c(n_in,n_in) ), c(n_in,2) ) )
vertexFCR <- rbind( vertexFC[1,]*cos(rot) + vertexFC[2,]*sin(rot) , vertexFC[1,]*sin(rot) + vertexFC[2,]*cos(rot) )

#save file centred
outTemp <- round( cbind( nodes, t( vertexFC ) ), 3 )
outFileName <- sprintf( '%s_axis_cent.dset', outName )
write.table( outTemp, row.name=FALSE, col.names=FALSE , outFileName )

#save file rotated
outTemp <- round( cbind( nodes, t( vertexFCR ) ), 3 )
outFileName <- sprintf( '%s_axis_cent_rot.dset', outName )
write.table( outTemp, row.name=FALSE, col.names=FALSE , outFileName )
