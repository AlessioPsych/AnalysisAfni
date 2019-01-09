##!/usr/bin/env bash

MPRAGE=$1
SEGMENTATIONFILE=$2
BOUNDARIESFILE=$3

hemisphereSeparation.sh $MPRAGE $SEGMENTATIONFILE /usr/local/bin/TT_icbm452+tlrc 1

3dTcat -prefix _del_boundariesSel.nii.gz $BOUNDARIESFILE[1,3,6] 

3dcalc -a _del_boundariesSel.nii.gz -b outputVolume.nii.gz -expr 'a*step(b)+not(step(b))*100' -prefix boundariesLeft.nii.gz
3dcalc -a _del_boundariesSel.nii.gz -b outputVolume.nii.gz -expr 'a*not(step(b))+step(b)*100' -prefix boundariesRight.nii.gz
rm _del_boundariesSel.nii.gz

defineBoundaries.sh boundariesLeft.nii.gz 1 2 1
mv boundariesThr.nii.gz boundariesThrLeft.nii.gz

defineBoundaries.sh boundariesRight.nii.gz 1 2 1
mv boundariesThr.nii.gz boundariesThrRight.nii.gz

generateSurfacesFromBoundaries.sh boundariesThrLeft.nii.gz 46 0-1-2 1500 1
mv surfaces_folder/ surfaces_folder_left/

generateSurfacesFromBoundaries.sh boundariesThrRight.nii.gz 46 0-1-2 1500 1
mv surfaces_folder/ surfaces_folder_right/

#instr <- 'binarizedConv.sh boundary01 surfaces_folder_left/'
#system( instr )
#system('mv convexityMap.dset convexityMap_left.dset')

#instr <- 'binarizedConv.sh boundary01 surfaces_folder_right/'
#system( instr )
#system('mv convexityMap.dset convexityMap_right.dset')


#boundariesFile <- 'V6690_GL_pdCorrectRegularT1_noBlur_stripped_alcortexMask_box_seg_boundaries.nii.gz[2,8,18]'

#instr <- sprintf('3dcalc -a %s -b outputVolume.nii.gz -expr \u0027a*step(b)+not(step(b))*100\u0027 -prefix boundariesLeft.nii.gz', boundariesFile)
# in matlab: instr <- sprintf('3dcalc -a %s -b outputVolume.nii.gz -expr ''a*step(b)+not(step(b))*100'' -prefix boundariesLeft.nii.gz', boundariesFile)

#print( instr )
#system( instr )

#instr <- sprintf('3dcalc -a %s -b outputVolume.nii.gz -expr \u0027a*not(step(b))+step(b)*100\u0027 -prefix boundariesRight.nii.gz', boundariesFile)
#print( instr )
#system( instr )

#instr <- sprintf('3dresample -dxyz 1 1 1 -prefix boundariesLeft_sm.nii.gz -input boundariesLeft.nii.gz -rmode Lin')
#print( instr )
#system( instr )

#instr <- sprintf('3dresample -dxyz 1 1 1 -prefix boundariesRight_sm.nii.gz -input boundariesRight.nii.gz -rmode Lin')
#print( instr )
#system( instr )

#instr <- sprintf('defineBoundaries.sh boundariesLeft_sm.nii.gz 1 2 1')
#print( instr )
#system( instr )
#system('mv boundariesThr.nii.gz boundariesThrLeft.nii.gz')

#instr <- '3dmask_tool -input boundariesThrLeft.nii.gz[0] -prefix boundariesThrLeft_dil01.nii.gz -dilate_input 1 -1'
#system( instr )
#instr <- '3dmask_tool -input boundariesThrLeft.nii.gz[1] -prefix boundariesThrLeft_dil02.nii.gz -dilate_input 1 -1'
#system( instr )
#instr <- '3dmask_tool -input boundariesThrLeft.nii.gz[2] -prefix boundariesThrLeft_dil03.nii.gz -dilate_input 1 -1'
#system( instr )
#instr <- '3dTcat -prefix boundariesThrLeft_dil.nii.gz boundariesThrLeft_dil01.nii.gz boundariesThrLeft_dil02.nii.gz boundariesThrLeft_dil03.nii.gz'
#system( instr )

#instr <- sprintf('defineBoundaries.sh boundariesRight_sm.nii.gz 1 2 1')
#print( instr )
#system( instr )
#system('mv boundariesThr.nii.gz boundariesThrRight.nii.gz')

#instr <- '3dmask_tool -input boundariesThrRight.nii.gz[0] -prefix boundariesThrRight_dil01.nii.gz -dilate_input 1 -1'
#system( instr )
#instr <- '3dmask_tool -input boundariesThrRight.nii.gz[1] -prefix boundariesThrRight_dil02.nii.gz -dilate_input 1 -1'
#system( instr )
#instr <- '3dmask_tool -input boundariesThrRight.nii.gz[2] -prefix boundariesThrRight_dil03.nii.gz -dilate_input 1 -1'
#system( instr )
#instr <- '3dTcat -prefix boundariesThrRight_dil.nii.gz boundariesThrRight_dil01.nii.gz boundariesThrRight_dil02.nii.gz boundariesThrRight_dil03.nii.gz'
#system( instr )

#instr <- sprintf('generateSurfacesFromBoundaries.sh boundariesThrLeft_dil.nii.gz 46 0-1-2 1500 1')
#system( instr )
#system('mv surfaces_folder/ surfaces_folder_left/')

#instr <- sprintf('generateSurfacesFromBoundaries.sh boundariesThrRight_dil.nii.gz 46 0-1-2 1500 1')
#system( instr )
#system('mv surfaces_folder/ surfaces_folder_right/')

#instr <- 'binarizedConv.sh boundary01 surfaces_folder_left/'
#system( instr )
#system('mv convexityMap.dset convexityMap_left.dset')

#instr <- 'binarizedConv.sh boundary01 surfaces_folder_right/'
#system( instr )
#system('mv convexityMap.dset convexityMap_right.dset')


#create ROI command
#'surf2vol.sh leftV1.1D.roi boundary00 V6690_GL_pdCorrectRegularT1_noBlur_stripped_alcortexMask_box_seg_depth.nii.gz 50 surfaces_folder_left/'

#visualize surface
# vglconnect compute-01
#'vglrun afniSurface.sh surfaces_folder_right/ V6690_GL_pdCorrectRegularT1_noBlur_stripped_alcortexMask_box_seg_anatomy.nii.gz'
