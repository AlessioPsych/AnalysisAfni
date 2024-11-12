


instr <- 'docker run -v /media/alessiofracasso/DATADRIVE1/testDocker/testMrtrix3:/testMrtrix3 --rm -it mrtrix3/mrtrix3 dwidenoise -nthreads 4 testMrtrix3/epiToTest_resampled.nii testMrtrix3/epiToTest_resampled_denoised.nii'

system( instr )




