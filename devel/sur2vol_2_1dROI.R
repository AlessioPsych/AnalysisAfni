
args <- commandArgs(T)
print( args )

#args <- c('leftCalcarine_surface_roi.dset','leftCalcarine_single_roi.1D')


surfaceFileName <- args[1]
roiName <- args[2]

surfaceTable <- read.table( surfaceFileName )

head( surfaceTable )

roiLineIdx <- which( surfaceTable[,7]==1 )

roiNodes <- surfaceTable[ roiLineIdx, 1 ]

roiTable <- cbind( roiNodes, rep( 1, length(roiNodes) ) )

write.table( roiTable, file = roiName, row.names = FALSE, col.names = FALSE )


