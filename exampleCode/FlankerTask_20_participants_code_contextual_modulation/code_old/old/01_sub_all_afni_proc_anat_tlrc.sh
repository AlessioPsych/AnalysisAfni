#!/usr/bin/env tcsh

#!/bin/tcsh

set foldername = "derivatives/processing_afni"

if ( -d $foldername ) then
    echo "Directory '$foldername' exists. Removing and recreating..."
    rm -rf $foldername
else
    echo "Directory '$foldername' does not exist. Creating..."
endif

mkdir $foldername

for subj in 'cat subjList_short.txt'; do

	# assign output directory name
	set output_dir = $subj.results
	mkdir $output_dir

done




