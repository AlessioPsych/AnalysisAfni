#!/bin/tcsh -xef

set atlasset = '/media/alessiofracasso/DATADRIVE1/Flanker/derivatives/processing_afni/sub-01.results/MNI152_2009_template.nii.gz'

set maskdset = '/media/alessiofracasso/DATADRIVE1/Flanker/derivatives/processing_afni/sub-01.results/mask_group+tlrc'

set dirA = '/media/alessiofracasso/DATADRIVE1/Flanker/derivatives/processing_afni'

set resultsdir = 'derivatives/test.results'

if ( -d $resultsdir ) then
    echo "Directory '$resultsdir' exists. Removing and recreating..."
    rm -rf $resultsdir
else
    echo "Directory '$resultsdir' does not exist. Creating..."
endif

mkdir $resultsdir

3dttest++ -prefix $resultsdir/Flanker_Inc_Cong_ttest \
	-mask $maskdset	                 \
	-setA Inc-Con			 \
	01 "$dirA/sub-01.results/stats.sub-01+tlrc[7]" \
	02 "$dirA/sub-02.results/stats.sub-02+tlrc[7]" \
	03 "$dirA/sub-03.results/stats.sub-03+tlrc[7]" \
	04 "$dirA/sub-04.results/stats.sub-04+tlrc[7]" \
	05 "$dirA/sub-05.results/stats.sub-05+tlrc[7]" \
	06 "$dirA/sub-06.results/stats.sub-06+tlrc[7]" \
	07 "$dirA/sub-07.results/stats.sub-07+tlrc[7]" \
	08 "$dirA/sub-08.results/stats.sub-08+tlrc[7]" \
	09 "$dirA/sub-09.results/stats.sub-09+tlrc[7]" \
	10 "$dirA/sub-10.results/stats.sub-10+tlrc[7]" \
	11 "$dirA/sub-11.results/stats.sub-11+tlrc[7]" \
	12 "$dirA/sub-12.results/stats.sub-12+tlrc[7]" \
	13 "$dirA/sub-13.results/stats.sub-13+tlrc[7]" \
	14 "$dirA/sub-14.results/stats.sub-14+tlrc[7]" \
	15 "$dirA/sub-15.results/stats.sub-15+tlrc[7]" \
	16 "$dirA/sub-16.results/stats.sub-16+tlrc[7]" \
	17 "$dirA/sub-17.results/stats.sub-17+tlrc[7]" \
	18 "$dirA/sub-18.results/stats.sub-18+tlrc[7]" \
	19 "$dirA/sub-19.results/stats.sub-19+tlrc[7]" \
	20 "$dirA/sub-20.results/stats.sub-20+tlrc[7]" \
	21 "$dirA/sub-21.results/stats.sub-21+tlrc[7]" \
	22 "$dirA/sub-22.results/stats.sub-22+tlrc[7]" \
	23 "$dirA/sub-23.results/stats.sub-23+tlrc[7]" \
	24 "$dirA/sub-24.results/stats.sub-24+tlrc[7]" \
	25 "$dirA/sub-25.results/stats.sub-25+tlrc[7]" \
	26 "$dirA/sub-26.results/stats.sub-26+tlrc[7]"


3dMEMA -prefix $resultsdir/Flanker_Inc_Cong_MEMA \
	-mask $maskdset	                 \
	-set Inc-Con			 \
	01 "$dirA/sub-01.results/stats.sub-01_REML+tlrc[7]" \
			"$dirA/sub-01.results/stats.sub-01_REML+tlrc[8]" \
	02 "$dirA/sub-02.results/stats.sub-02_REML+tlrc[7]" \
			"$dirA/sub-02.results/stats.sub-02_REML+tlrc[8]" \
	03 "$dirA/sub-03.results/stats.sub-03_REML+tlrc[7]" \
			"$dirA/sub-03.results/stats.sub-03_REML+tlrc[8]" \
	04 "$dirA/sub-04.results/stats.sub-04_REML+tlrc[7]" \
			"$dirA/sub-04.results/stats.sub-04_REML+tlrc[8]" \
	05 "$dirA/sub-05.results/stats.sub-05_REML+tlrc[7]" \
			"$dirA/sub-05.results/stats.sub-05_REML+tlrc[8]" \
	06 "$dirA/sub-06.results/stats.sub-06_REML+tlrc[7]" \
			"$dirA/sub-06.results/stats.sub-06_REML+tlrc[8]" \
	07 "$dirA/sub-07.results/stats.sub-07_REML+tlrc[7]" \
			"$dirA/sub-07.results/stats.sub-07_REML+tlrc[8]" \
	08 "$dirA/sub-08.results/stats.sub-08_REML+tlrc[7]" \
			"$dirA/sub-08.results/stats.sub-08_REML+tlrc[8]" \
	09 "$dirA/sub-09.results/stats.sub-09_REML+tlrc[7]" \
			"$dirA/sub-09.results/stats.sub-09_REML+tlrc[8]" \
	10 "$dirA/sub-10.results/stats.sub-10_REML+tlrc[7]" \
			"$dirA/sub-10.results/stats.sub-10_REML+tlrc[8]" \
	11 "$dirA/sub-11.results/stats.sub-11_REML+tlrc[7]" \
			"$dirA/sub-11.results/stats.sub-11_REML+tlrc[8]" \
	12 "$dirA/sub-12.results/stats.sub-12_REML+tlrc[7]" \
			"$dirA/sub-12.results/stats.sub-12_REML+tlrc[8]" \
	13 "$dirA/sub-13.results/stats.sub-13_REML+tlrc[7]" \
			"$dirA/sub-13.results/stats.sub-13_REML+tlrc[8]" \
	14 "$dirA/sub-14.results/stats.sub-14_REML+tlrc[7]" \
			"$dirA/sub-14.results/stats.sub-14_REML+tlrc[8]" \
	15 "$dirA/sub-15.results/stats.sub-15_REML+tlrc[7]" \
			"$dirA/sub-15.results/stats.sub-15_REML+tlrc[8]" \
	16 "$dirA/sub-16.results/stats.sub-16_REML+tlrc[7]" \
			"$dirA/sub-16.results/stats.sub-16_REML+tlrc[8]" \
	17 "$dirA/sub-17.results/stats.sub-17_REML+tlrc[7]" \
			"$dirA/sub-17.results/stats.sub-17_REML+tlrc[8]" \
	18 "$dirA/sub-18.results/stats.sub-18_REML+tlrc[7]" \
			"$dirA/sub-18.results/stats.sub-18_REML+tlrc[8]" \
	19 "$dirA/sub-19.results/stats.sub-19_REML+tlrc[7]" \
			"$dirA/sub-19.results/stats.sub-19_REML+tlrc[8]" \
	20 "$dirA/sub-20.results/stats.sub-20_REML+tlrc[7]" \
			"$dirA/sub-20.results/stats.sub-20_REML+tlrc[8]" \
	21 "$dirA/sub-21.results/stats.sub-21_REML+tlrc[7]" \
			"$dirA/sub-21.results/stats.sub-21_REML+tlrc[8]" \
	22 "$dirA/sub-22.results/stats.sub-22_REML+tlrc[7]" \
			"$dirA/sub-22.results/stats.sub-22_REML+tlrc[8]" \
	23 "$dirA/sub-23.results/stats.sub-23_REML+tlrc[7]" \
			"$dirA/sub-23.results/stats.sub-23_REML+tlrc[8]" \
	24 "$dirA/sub-24.results/stats.sub-24_REML+tlrc[7]" \
			"$dirA/sub-24.results/stats.sub-24_REML+tlrc[8]" \
	25 "$dirA/sub-25.results/stats.sub-25_REML+tlrc[7]" \
			"$dirA/sub-25.results/stats.sub-25_REML+tlrc[8]" \
	26 "$dirA/sub-26.results/stats.sub-26_REML+tlrc[7]" \
			"$dirA/sub-26.results/stats.sub-26_REML+tlrc[8]"


cp $atlasset $resultsdir 
