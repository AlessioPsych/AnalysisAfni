# to run ' tcsh 04_00_forLoop_checkProcessing_ORIG.sh '

#!/bin/tcsh -xef

set maindir = "/mnt/disk01/ds000102_FlankerTask"
set codedir = ${maindir}/code

cd $maindir

# creates output folder
set output_folder = "derivatives/images_afni_ORIG_from_denoised_data"
set input_folder = "derivatives/processing_afni_denoised_incongruent_congruent"

if ( -d $output_folder ) then
    echo "Directory '$output_folder' exists. Removing and recreating..."
    rm -rf $output_folder
else
    echo "Directory '$output_folder' does not exist. Creating..."
endif

mkdir $output_folder

# creates subjList.txt file
if (! -f subjList.txt) then
    ls | grep ^sub- > subjList.txt
endif


# runs code across participants
cd $maindir
foreach i (`cat subjList.txt`)	
    cd $maindir
    cd $codedir
    tcsh 04_01_visualizeProcessed_ORIG.sh $i $output_folder $input_folder $maindir
end




