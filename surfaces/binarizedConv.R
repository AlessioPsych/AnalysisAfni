args <- commandArgs(T)
print( args )



system( sprintf( 'SurfaceMetrics -conv -i_1D %s/%s_sm.1D.coord %s/%s_or.1D.topo -prefix conv_data', args[2], args[1], args[2], args[1] ) )
system( sprintf('SurfSmooth -spec %s/spec.surfaces.smoothed -surf_A %s/%s_sm.1D.coord  -met HEAT_05 -input conv_data.conv.1D.dset -fwhm 15 -output conv_data.conv_smooth.1D.dset', args[2], args[2], args[1], args[1] ) )
convMap <- read.table('conv_data.conv_smooth.1D.dset', as.is=TRUE )
convMapThr01 <- ifelse( convMap[,2] < 0, 0.8, 0.2 )
convMapThr02 <- ifelse( convMap[,2] > 0, 0.8, 0.2 )
convMapBin <- cbind( convMap, convMapThr01, convMapThr02  )
write.table( convMapBin, file='convexityMap.dset', row.names=FALSE, col.names=FALSE)
system('rm conv_data.conv_smooth.1D.dset')
system('rm conv_data.conv.niml.dset')
system('rm conv_data.conv.1D.dset')
