
#using docker and R (or matlab)
#start terminal with admin privileges:

sudo -s

#start R

R

instr <- 'docker run -v /media/alessiofracasso/DATADRIVE1/testDocker/testMrtrix3:/testMrtrix3 --rm -it mrtrix3/mrtrix3'
system( instr )
cd testMrtrix3/
dwidenoise epiToTest.nii epiDenoisedTestR.nii
chmod 777 epiDenoisedTestR.nii
exit



