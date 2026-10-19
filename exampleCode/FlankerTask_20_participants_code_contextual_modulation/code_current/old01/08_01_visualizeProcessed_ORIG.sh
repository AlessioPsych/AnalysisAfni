#!/bin/tcsh -xef

set subj  = $argv[1]
#set subj  = 'sub-01'
set maindir = $argv[4]
set targetdir = $argv[2]
set datadir = $argv[3]

echo $subj
echo $maindir
echo $targetdir
echo $datadir

cd $maindir
cd $datadir
cd {$subj}.results/

@chauffeur_afni                     \
        -ulay    anat_final.{$subj}+orig         \
        -prefix  ${maindir}/${targetdir}/{$subj}-anatOrig         \
        -montx 3 -monty 5               \
        -set_xhairs OFF                 \
        -label_mode 0 -label_size 0     \
        -do_clean  


3dTcat -prefix delMeTemp.nii.gz pb02.{$subj}.r01.volreg+orig'[1]'

@chauffeur_afni                     \
        -ulay    delMeTemp.nii.gz         \
        -prefix  ${maindir}/${targetdir}/{$subj}-volregOrig         \
        -montx 3 -monty 5               \
        -set_xhairs OFF                 \
        -label_mode 0 -label_size 0     \
        -do_clean  

rm delMeTemp.nii.gz


