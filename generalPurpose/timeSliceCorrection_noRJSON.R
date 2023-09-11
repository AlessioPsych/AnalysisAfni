args <- commandArgs(T)
print( args )

#to debug
#setwd('/media/alessiof/Data/tests/SEF_visual_response/HMR28')
#args <- c('EPI/','*.nii','EPI_jsons/','EPI_tsc/')

mainDir <- getwd()
generalPurposeDir <- Sys.getenv( x='AFNI_TOOLBOXDIRGENERALPURPOSE' )
afniInstallDir <- Sys.getenv( x='AFNI_INSTALLDIR' ) 


#arrange input
inputEpiDir <- args[1]
epiFormat <- args[2]
inputJsonDir <- args[3]
outputNameDir <- args[4]

print( sprintf('input EPI directory = %s', inputEpiDir ) )
print( sprintf('EPI file format = %s', epiFormat ) )
print( sprintf('input JSON directory = %s', inputJsonDir ) )
print( sprintf('output directory = %s', outputNameDir ) )

# create output directory
instr <- sprintf('mkdir %s', outputNameDir ); system( instr )

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


# get epi filenames
setwd( inputEpiDir )
filesWorkDir <- dir(pattern=epiFormat);
setwd( mainDir )

# get json filenames
setwd( inputJsonDir )
filesJsons <- dir(pattern='*.json');
setwd( mainDir )

# check that wehave the same number of epi and json files, abort otherwise
flagCheckNFiles <- length( filesJsons ) == length( filesWorkDir )
stopifnot( flagCheckNFiles )

# perform ts correction
for (k in 1:length(filesWorkDir)) {
  
  print( sprintf('processing epi file: %s', filesWorkDir[k] ) )

  print( sprintf('loading json file: %s', filesJsons[k] ) )  
  # Specify the path to your JSON file
  json_file <- sprintf( '%s%s', inputJsonDir, filesJsons[k] )

  # Read the JSON file as text lines
  json_lines <- readLines(json_file)
	
  # Parse the JSON content
  parsed_json <- parse_json(json_lines)

  # Access the first element of the JSON array
  first_element <- parsed_json[[1]]
       
  print( sprintf('writing slice timing for file: %s', filesJsons[k] ) )
  sliceTimingData <- first_element$SliceTiming
  
  write.table( x=t(sliceTimingData), file = sprintf( '%s_ttt_slice_timing_file.1D', inputEpiDir ), row.names = FALSE, col.names = FALSE  )
  
  instr <- sprintf('3dTshift -tpattern @%s_ttt_slice_timing_file.1D -prefix %s%s %s%s', inputEpiDir, outputNameDir, filesWorkDir[k], inputEpiDir, filesWorkDir[k]  ); system(instr)
  
  # clean up
  instr <- sprintf( 'rm %s_ttt_slice_timing_file.1D', inputEpiDir ); system( instr )
  
}
