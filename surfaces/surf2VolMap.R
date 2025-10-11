args <- commandArgs(T)
print( args )

source( sprintf('%s/AFNIio.R', args[4] ) )
source( sprintf('%s/coordinateFromLinearIndex.r', args[5] ) )
source( sprintf('%s/linearIndexFromCoordinate.r', args[5] ) )
library(RANN)

commandLine <- sprintf( '3dSurf2Vol -spec %sspec.surfaces.smoothed -surf_A %s%s_sm.1D.coord -sv %s -grid_parent %s -map_func %s -prefix %s -sdata_1D %s', args[7], args[7], args[2], args[3], args[3], args[8], 'surfVolMap', args[1] )

print( commandLine )
system( commandLine )


