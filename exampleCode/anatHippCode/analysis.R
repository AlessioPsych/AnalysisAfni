mainDir <- 'C:/Users/alessiof/Downloads'
setwd( mainDir )
print( getwd() )
foldersToProcess <- dir('anatHippData')
for (n in 1:length(foldersToProcess)) { # n <- 1
  setwd( mainDir ) 
  print( getwd() )
  currentAnatomy <- foldersToProcess[n]
  instr <- paste( 'docker', 
  'run', 
  '-it', 
  '-v', 
  sprintf('%s/anatHippData/%s:/bids', mainDir, currentAnatomy), 
  '-v', 
  sprintf('%s/anatHippData/%s_output:/output', mainDir, currentAnatomy), 
  'khanlab/hippunfold:latest', 
  '/bids', 
  '/output participant', 
  '--modality T1w',
  '--cores 8' ) # -n --cores 8
  print( instr )
  system( instr )
}