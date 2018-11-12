#!/usr/bin/env bash

coordFiles=(*_sm.1D.coord)
topoFiles=(*.topo)

for (( i=0; i<=9; i++ )) do
        outputFile=$(printf 'surfMapping0%s.surf' $i)
	instr1=$(printf 'SurfToSurf -i_1D %s %s -i_1D %s %s -prefix %s' ${coordFiles[i]} ${topoFiles[i]} ${coordFiles[i+1]} ${topoFiles[i+1]} ${outputFile} )
        echo ${i}
	${instr1}	
done


#SurfToSurf -i_1D boundary08_sm.1D.coord boundary08.1D.topo \
#	-i_1D boundary09_sm.1D.coord boundary09.1D.topo \
#	-prefix thickness.surf

