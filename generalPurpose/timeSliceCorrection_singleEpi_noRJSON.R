args <- commandArgs(T)
print( args )

#to debug
#setwd('/media/alessiof/Data/tests/SEF_visual_response/HMR28')
#args <- c('EPI/','*.nii','EPI_jsons/','EPI_tsc/')

mainDir <- getwd()
#generalPurposeDir <- Sys.getenv( x='AFNI_TOOLBOXDIRGENERALPURPOSE' )
#afniInstallDir <- Sys.getenv( x='AFNI_INSTALLDIR' ) 

#arrange input
inputEpiFile <- args[1]
inputJsonFile <- args[2]
outputEpiFile <- args[3]

print( sprintf('input EPI file = %s', inputEpiFile) )
print( sprintf('input JSON file = %s', inputJsonFile ) )
print( sprintf('output File = %s', outputEpiFile ) )

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

print( sprintf('processing epi file: %s', inputEpiFile ) )

print( sprintf('loading json file: %s', inputJsonFile ) )  

# Read the JSON file as text lines
json_lines <- readLines(inputJsonFile)
	
# Parse the JSON content
parsed_json <- parse_json(json_lines)

# Access the first element of the JSON array
first_element <- parsed_json[[1]]
       
print( sprintf('writing slice timing for file: %s', inputJsonFile ) )
sliceTimingData <- first_element$SliceTiming
  
write.table( x=t(sliceTimingData), file = 'ttt_slice_timing_file.1D', row.names = FALSE, col.names = FALSE  )
  
instr <- sprintf('3dTshift -tpattern @ttt_slice_timing_file.1D -prefix %s %s', outputEpiFile, inputEpiFile ); 
system(instr)
  
# clean up
instr <- sprintf( 'rm ttt_slice_timing_file.1D' ); 
system( instr )
  

