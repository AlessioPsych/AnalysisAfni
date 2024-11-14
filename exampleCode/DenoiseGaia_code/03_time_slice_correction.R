rm( list=ls() ); gc();

#need to load module R, afni, and activate conda environment beforehand

#### Move to the desired folder ####
sprintf('Move to the desired folder:')
mainDir <- '/scratch/af4887/Proj_Gaia_David'
inputDir <- '/scratch/af4887/Proj_Gaia_David/rawdata_denoised'
outputDir <- '/scratch/af4887/Proj_Gaia_David/EPI_timeSlicedCorrected'
setwd( mainDir )
sprintf( 'Current folder: %s', print( getwd() ) )
singleSubjectFolders <- dir( inputDir )
runCodeFlag <- 1

# clean up output folder, if it exists
if ( dir.exists( outputDir ) ) {
  instr <- sprintf('rm -R %s', outputDir )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  instr <- sprintf('mkdir %s', outputDir )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
}else{
  instr <- sprintf('mkdir %s', outputDir )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
}

# function to parse json file
parse_json <- function(json_lines) {
  # Initialize an empty list to store the JSON objects
  json_list <- list()
  
  # Initialize an empty character string to store the current JSON object
  current_object <- NULL
  
  # Loop through each line of JSON content
  for (line in json_lines) {
    # Append the current line to the current object string
    current_object <- paste(current_object, line, sep = "")
    
    # Check if the current object is complete (ends with '}' character)
    if (grepl("\\}$", current_object)) {
      # Parse the JSON object and add it to the list
      json_list <- c(json_list, list(jsonlite::fromJSON(current_object)))
      
      # Reset the current object
      current_object <- NULL
    }
  }
  
  return(json_list)
}


for ( nSubj in 1 : length( singleSubjectFolders )  ) {# nSubj <- 1 length( singleSubjectFolders )
  
  sprintf('Move to the individual participant folder:')
  setwd( inputDir )
  setwd( singleSubjectFolders[ nSubj ] )  
  sprintf( 'Current folder: %s', print( getwd() ) )
  setwd( 'func' )
  json_files_toProcess <- dir( pattern='*bold.json' )
  nifti_files_toProcess <- dir( patter='*bold.nii*' )
  
  # check that wehave the same number of epi and json files, abort otherwise
  flagCheckNFiles <- length( json_files_toProcess ) == length( nifti_files_toProcess )
  stopifnot( flagCheckNFiles )
  
  # make subject specific output folder
  instr <- sprintf('mkdir %s/%s', outputDir, singleSubjectFolders[ nSubj ] )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  
  for ( nJson in 1:length( json_files_toProcess ) ) { #nJson <- 1 length( json_files_toProcess )
  
    print( sprintf('loading json file: %s', json_files_toProcess[ nJson ] ) )  
    print( sprintf('to process with nifti file: %s', nifti_files_toProcess[ nJson ] ) ) 
    outPutFileName <- paste( strsplit( nifti_files_toProcess[ nJson ], '[.]' )[[1]][1], '_TScorr.nii.gz', sep='' )
    json_file <- json_files_toProcess[ nJson ]
    
    # Read the JSON file as text lines
    json_lines <- suppressWarnings( readLines(json_file) )
    
    # Parse the JSON content
    parsed_json <- parse_json(json_lines)
    
    # Access the first element of the JSON array
    first_element <- parsed_json[[1]]
    
    sliceTimingData <- first_element$SliceTiming
    
    print( sprintf('writing slice timing for file: %s', json_files_toProcess[ nJson ] ) )
    if (runCodeFlag==1) { 
      write.table( x=t(sliceTimingData), file = sprintf( '%s/%s/_ttt_slice_timing_file.1D', outputDir, singleSubjectFolders[ nSubj ] ), row.names = FALSE, col.names = FALSE  )
    }

    print( sprintf('perform slice timing for file: %s', nifti_files_toProcess[ nJson ] ) )
    instr <- sprintf('afni 3dTshift -tpattern @%s/%s/_ttt_slice_timing_file.1D -prefix %s/%s/%s %s', 
                     outputDir, 
                     singleSubjectFolders[ nSubj ], 
                     outputDir, 
                     singleSubjectFolders[ nSubj ],
                     outPutFileName,
                     nifti_files_toProcess[ nJson ] )
    print( instr )
    if (runCodeFlag==1) { 
      system(instr)
    }  
    
    # clean up
    instr <- sprintf( 'rm %s/%s/_ttt_slice_timing_file.1D', outputDir, singleSubjectFolders[ nSubj ] )
    print( instr )
    if (runCodeFlag==1) { 
      system(instr)
    }  
    
  }
  
}

setwd( mainDir )
sprintf( 'done, current folder: %s', print( getwd() ) )

