#!/bin/tcsh -xef

set maindir = "/media/alessiofracasso/DATADRIVE1/Flanker/"
set targetdir = "derivatives/Freesurfer_output"

echo $maindir
echo $targetdir

# creates subjList.txt file
if (! -f subjList.txt) then
    ls | grep ^sub- > subjList.txt
endif

# runs code across participants
foreach i (`cat subjListShort.txt`)

	set subj = $i
	cd $maindir
	cd $targetdir
	cd $subj
	cd Freesurfer_result
	cd SUMA

	echo $PWD
	echo $subj
	
	# combine left and right upsampled ribbons
	if ( -f gray_matter_ribbon.nii.gz ) then 
		rm gray_matter_ribbon.nii.gz 
	endif	
	if (! -f gray_matter_ribbon.nii.gz) then	
		3dcalc -a rh.ribbon.nii -b lh.ribbon.nii -expr 'a+b' -prefix gray_matter_ribbon.nii.gz
	endif

	# get white matter mask from aparc+aseg segmentation
	if ( -f white_matter_mask.nii.gz ) then 
		rm white_matter_mask.nii.gz 
	endif
	if (! -f white_matter_mask.nii.gz) then
		3dcalc -a aparc+aseg.nii -expr 'within(a,1.5,2.5)+within(a,40.5,41.5)' -prefix white_matter_mask.nii.gz
	endif
	
	# combine white matter mask and gray matter ribbon at original resolution
	if ( -f gray_matter_mask_out.nii.gz ) then 
		rm gray_matter_mask_out.nii.gz 
	endif
	if (! -f gray_matter_mask_out.nii.gz) then
		3dcalc -a white_matter_mask.nii.gz -b gray_matter_ribbon.nii.gz -expr 'a+b' -prefix gray_matter_mask_out.nii.gz
	endif
end


