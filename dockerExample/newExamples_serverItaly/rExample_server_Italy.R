


instr <- 'docker run -v /home/fracasso/Data/test/testMrtrixDenoise:/testMrtrix3 --rm -it mrtrix3/mrtrix3 dwidenoise -nthreads 8 testMrtrix3/07_testNifti.nii testMrtrix3/07_testNifti_denoised_github.nii'

system( instr )




