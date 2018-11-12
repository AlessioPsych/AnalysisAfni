#!/usr/bin/env bash

coordFiles=(boundary*_sm.1D.coord)
topoFiles=(boundary*.topo)

    for (( i=5; i<${#coordFiles[@]}-1; i++ )) do
    outputFile=$(printf 'surfMapping%s_01.surf' $i)
    part=$(printf 'SurfToSurf -i_1D %s %s -i_1D %s %s -prefix %s' ${topoFiles[i]} ${coordFiles[i]} ${topoFiles[i+1]} ${coordFiles[i+1]} ${outputFile} )
    echo ${i}
    echo ${part}
    #${part}
    #echo ${outputFile}
done

for (( i=5; i>0; i-- )) do
    outputFile=$(printf 'surfMapping%s_02.surf' $i)
    part=$(printf 'SurfToSurf -i_1D %s %s -i_1D %s %s -prefix %s' ${topoFiles[i]} ${coordFiles[i]} ${topoFiles[i-1]} ${coordFiles[i-1]} ${outputFile} )
    echo ${i}
    echo ${part}
    #${part}
    #echo ${outputFile}
done

