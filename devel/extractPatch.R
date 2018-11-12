args <- commandArgs(T)
print( args )

args <- c( 'boundary03', 'test08.1D.roi', '0' )

boundary <- args[1]
roiFile <- args[2]
roiColumn <- args[3]

outputName <- strsplit( roiFile, '[.]' )[[1]][1]
outputName <- sprintf( '%s_patch', outputName )
instr  <- sprintf( 'SurfPatch -spec surfaces_folder/spec.surfaces.smoothed -surf surfaces_folder/%s_sm.1D.coord -input %s %s -1 -prefix %s -patch2surf -adjust_contour -hits 3 -out_type 1D', boundary, roiFile,  roiColumn, outputName )
system( instr )