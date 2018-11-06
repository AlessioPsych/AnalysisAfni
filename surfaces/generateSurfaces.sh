#!/usr/bin/env bash

if [ -z "$1" ]
 then
  echo 'Inputs:'
  echo 'NSURFACES=$1'
  echo 'SMOOTHING=$2'
  echo 'INFLATEITER=$3'
  echo 'INFLATEINDEX=$4' 
  exit 1
fi


NSURFACES=$1
SMOOTHING=$2
INFLATEITER=$3
INFLATEINDEX=$4



## generate all isosurfaces (unsmoothed)
for (( i=0; i<=($NSURFACES-1); i++ ))
do
    outputFileName=$(printf 'boundariesThr.nii.gz[%s]' $i)
    if [ "$i" -lt 10 ] 
    then
        instr1=$(printf 'IsoSurface -input %s -isoval 1 -o_vec boundary0%s_or' $outputFileName $i)
    else
        instr1=$(printf 'IsoSurface -input %s -isoval 1 -o_vec boundary%s_or' $outputFileName $i)
    fi
        echo $i
        echo ${outputFileName}
        ${instr1}
done

## prepare spec file (unsmoothed)
coordFiles=(*.coord)
topoFiles=(*.topo)
instr2='quickspec -spec spec.surfaces'
for (( i=0; i<=($NSURFACES-1); i++ )) do
    part=$(printf ' -tn 1D %s %s '${topoFiles[i]} ${coordFiles[i]} )
    instr2=$(printf ' %s %s ' ${instr2} ${part} )
    echo ${i}
    echo ${instr2}
done

${instr2}

## smooth surfaces
for (( i=0; i<=($NSURFACES-1); i++ ))
do
    if (( $i < 10 )); then
        instr3=$(printf 'SurfSmooth -spec spec.surfaces -surf_A %s -met LM -Niter %s -output boundary0%s_sm.1D.coord' ${coordFiles[i]} $SMOOTHING $i )
    else
        instr3=$(printf 'SurfSmooth -spec spec.surfaces -surf_A %s -met LM -Niter %s -output boundary%s_sm.1D.coord' ${coordFiles[i]} $SMOOTHING $i )
    fi
    echo ${i}
    echo ${instr3}
    $instr3 
done

## inflate surface
selectedSurface=${coordFiles[$INFLATEINDEX]}
SurfSmooth  -spec spec.surfaces -surf_A $selectedSurface -met NN_geom \
                   -output 'inflated_sm.1D.coord' -Niter $INFLATEITER \
                   -match_area 0.01
cp ${topoFiles[$INFLATEINDEX]} 'inflated.1D.topo'
#SurfSmooth  -spec spec.surfaces -surf_A $selectedSurface -met LM \
#                   -output 'inflated_sm.1D.coord' -Niter $INFLATEITER
#cp ${topoFiles[0]} 'inflated.1D.topo'

## prepare spec file (smoothed)
coordFiles=(*_sm.1D.coord)
topoFiles=(*.topo)
instr4='quickspec -spec spec.surfaces.smoothed'
for (( i=0; i<=($NSURFACES); i++ )) 
do
	if (( $i < $NSURFACES )); then	
		    part=$(printf ' -tn 1D %s %s '${topoFiles[i]} ${coordFiles[i]} )
		    instr4=$(printf ' %s %s ' ${instr4} ${part} )
	else
		    part=$(printf ' -tsnad 1D S_inflated %s %s N %s '${coordFiles[$INFLATEINDEX]} ${coordFiles[i]} ${topoFiles[i]}  )
		    instr4=$(printf ' %s %s ' ${instr4} ${part} )
	fi
	echo ${i}
	echo ${instr4}
done

${instr4} 

rm *_or.1D.coord
rm spec.surfaces

echo 'INFLATED SURFACE:'
echo $selectedSurface

mkdir surfaces_folder
mv *.1D.topo surfaces_folder/
mv *.1D.coord surfaces_folder/
mv spec.surfaces.smoothed surfaces_folder/


## clean up
#rm *.coord
#rm *.topo
#rm surfaces
#rm *.1D
#rm boundariesThr.nii.gz
#rm spec.surfaces.smoothed


#afni -niml &
#suma -spec spec.surfaces.smoothed -sv occipital_pole.nii.gz &




#IsoSurface -input GM_CSF_boundary.nii.gz -isoval 1 -o_vec CSF

#quickspec -tn 1D WM.1D.coord WM.1D.topo \
#	  -tn 1D CSF.1D.coord CSF.1D.topo \
#	  -spec quick_surfaces 

#SurfSmooth -spec quick_surfaces  -surf_A WM.1D.coord -met LM -Niter 50 -output WM_sm.1D.coord

#SurfSmooth -spec quick_surfaces  -surf_A CSF.1D.coord -met LM -Niter 50 -output CSF_sm.1D.coord

#rm quick_surfaces 
#rm WM.1D.coord
#rm CSF.1D.coord

#quickspec -tn 1D WM_sm.1D.coord WM.1D.topo \
#	  -tn 1D CSF_sm.1D.coord CSF.1D.topo \
#	  -spec surfaces 

#3dVol2Surf                                    \
#       -spec         spec.surfaces.smoothed   \
#       -surf_A       boundary05_sm.1D.coord   \
#       -sv           occipital_pole.nii.gz    \
#       -grid_parent  occipital_pole.nii.gz    \
#       -map_func     mask                     \
#       -out_1D       boundary05_vals.1D

#3dVol2Surf                                    \
#       -spec         spec.surfaces.smoothed   \
#       -surf_A       boundary08_sm.1D.coord   \
#       -surf_B       boundary09_sm.1D.coord   \
#       -sv           occipital_pole.nii.gz    \
#       -grid_parent  occipital_pole.nii.gz    \
#      -map_func     ave                      \
#       -out_1D       WM_vals_1.1D

#SurfToSurf -i_1D boundary08_sm.1D.coord boundary08.1D.topo \
#	-i_1D boundary09_sm.1D.coord boundary09.1D.topo \
#	-prefix thickness.surf



#afni -niml &
#suma -spec surfaces -sv occipital_pole.nii.gz
