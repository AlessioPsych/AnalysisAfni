#!/bin/tcsh -xef

# creates output folder
set maindir="/home/fracasso/Data/openNeuro/ds000102"
set codedir="/code"
set foldername = "derivatives/images_processing_Freesurfer"
cd $maindir

if ( -d $foldername ) then
    echo "Directory '$foldername' exists. Removing and recreating..."
    rm -rf $foldername
else
    echo "Directory '$foldername' does not exist. Creating..."
endif

mkdir $foldername

# creates subjList.txt file
if (! -f subjList.txt) then
    ls | grep ^sub- > subjList.txt
endif

# runs code across participants
foreach i (`cat subjList.txt`)
    tcsh $maindir/$codedir/01_02_01_visualizeSegmentations.sh $i
end




