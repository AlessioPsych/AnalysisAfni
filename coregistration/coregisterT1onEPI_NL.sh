MPRAGE=$1
EPI=$2
MODE=$3
UNIFEPI=$4
UNIFMPRAGE=$5

3dAutomask -apply_prefix tt_epi+orig $EPI

if (( $UNIFEPI == 1 ))
	then
	3dUnifize -prefix tt_epi_unif+orig -T2 -GM tt_epi+orig
fi
if (( $UNIFEPI == 0 ))
	then
	3dcopy tt_epi+orig tt_epi_unif+orig
fi


if (( $UNIFMPRAGE == 1 ))
	then
	3dUnifize -prefix tt_MPRAGE_unif+orig -T2 -T2 -GM $MPRAGE
fi
if (( $UNIFMPRAGE == 2 ))
	then
	3dUnifize -prefix tt_MPRAGE_unif+orig -T2 -GM $MPRAGE
fi
if (( $UNIFMPRAGE == 3 ))
	then
	3dUnifize -prefix tt_MPRAGE_unif+orig $MPRAGE
fi
if (( $UNIFMPRAGE == 4 ))
	then
	3dcopy $MPRAGE tt_MPRAGE_unif+orig
fi


if (( $MODE == 1 )) 
	then
		3dQwarp -source tt_MPRAGE_unif+orig   \
		              -base tt_epi_unif+orig       \
		              -blur 0 0 \
			      -patchmin 19 \
			      -iwarp
	rm tt_*
	mv Qwarp+orig.BRIK Qwarp_1+orig.BRIK
	mv Qwarp+orig.HEAD Qwarp_1+orig.HEAD
	mv Qwarp_WARP+orig.BRIK Qwarp_WARP_1+orig.BRIK
	mv Qwarp_WARP+orig.HEAD Qwarp_WARP_1+orig.HEAD
	mv Qwarp_WARPINV+orig.BRIK Qwarp_WARPINV_1+orig.BRIK
	mv Qwarp_WARPINV+orig.HEAD Qwarp_WARPINV_1+orig.HEAD
fi

if (( $MODE == 2 )) 
	then
		3dQwarp -source tt_MPRAGE_unif+orig   \
              -base tt_epi_unif+orig       \
              -blur 0 0 \
              -plusminus \
	      -iwarp
	rm tt_*

fi

if (( $MODE == 3 )) 
	then
		3dQwarp -source tt_MPRAGE_unif+orig   \
		       -base tt_epi_unif+orig \
		       -blur 0 0 \
                       -iwarp
	rm tt_*
	mv Qwarp+orig.BRIK Qwarp_3+orig.BRIK
	mv Qwarp+orig.HEAD Qwarp_3+orig.HEAD
	mv Qwarp_WARP+orig.BRIK Qwarp_WARP_3+orig.BRIK
	mv Qwarp_WARP+orig.HEAD Qwarp_WARP_3+orig.HEAD
	mv Qwarp_WARPINV+orig.BRIK Qwarp_WARPINV_3+orig.BRIK
	mv Qwarp_WARPINV+orig.HEAD Qwarp_WARPINV_3+orig.HEAD
fi




