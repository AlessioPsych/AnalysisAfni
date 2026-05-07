# runs the mrtrix3 docker passing the '/media/alessiofracasso/DATADRIVE1/testDocker/testMrtrix3' into the folder '/testMrtrix3'
# it is possible to operate on the files in the '/testMrtrix3', then remember to free up the permissions of the newly created files, see line 7 of the current file

sudo docker run -v /media/alessiofracasso/DATADRIVE1/testDocker/testMrtrix3:/testMrtrix3 --rm -it mrtrix3/mrtrix3
cd testMrtrix3/
dwidenoise epiToTest.nii epiDenoised.nii
chmod 777 epiDenoised.nii

# nighres
# this line runs the docker in the command line, still working from here
sudo docker run --rm -it nighres
