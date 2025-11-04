#!/bin/tcsh -xef

set maindir = $PWD

# creates output folder
set foldername = "derivatives/images_processing_afni"

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
foreach i (`cat subjList_short.txt`)
	cd $maindir
	cd $foldername
	echo $PWD	
end

