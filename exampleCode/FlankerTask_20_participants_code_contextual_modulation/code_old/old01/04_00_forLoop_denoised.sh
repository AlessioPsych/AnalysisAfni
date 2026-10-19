#!/bin/tcsh -xef

set main_folder = '/home/fracasso/Data/openNeuro/ds000102'
set code_folder = $main_folder/code
set input_folder = $main_folder/derivatives/mrtrix3
set output_folder = $main_folder/derivatives/processing_afni_denoised
echo $main_folder
echo $code_folder
echo $input_folder
echo $output_folder

cd $main_folder

# creates output folder

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
foreach i (`cat subjList.txt`)
    echo instr: tcsh $code_folder/04_01_sub_xx_afni_proc_denoised.sh $i $main_folder $input_folder
    echo instr mv ${i}.results $output_folder
    tcsh $code_folder/04_01_sub_xx_afni_proc_denoised.sh $i $main_folder $input_folder
    mv ${i}.results $output_folder
    mv proc.${i} $output_folder
    mv output.proc.${i} $output_folder    
    echo ................
    echo ................
    echo ................
end




