function matCall(afniMatDir,myFile,myDir)

sprintf('%s',afniMatDir)
sprintf('%s',myFile)
sprintf('%s',myDir)
addpath(genpath(afniMatDir))

cd myDir

%system( sprintf( '3dcopy %s %s', myFile, 'p+orig' ) )

exit
