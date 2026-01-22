maindir="/home/fracasso/Data/openNeuro/ds000102/derivatives/mrtrix3"

echo "$maindir"

cd "$maindir"
echo "Current folder: $PWD"

#Check whether the file subjList.txt exists; if not, create it
if [ -f subjList.txt ]; then
	rm subjList.txt
fi
if [ ! -f subjList.txt ]; then
	ls | grep sub- > subjList.txt
fi

#Loop over all subjects and format timing files into FSL format
for subj in `cat subjList.txt`; do
	cd $subj/func

	cat ${subj}_task-flanker_run-1_events.tsv | awk '{if ($2=="2.0") {print $1, $4, 1}}' > RT_correctWrong_run1.txt

	cat ${subj}_task-flanker_run-2_events.tsv | awk '{if ($2=="2.0") {print $1, $4, 1}}' > RT_correctWrong_run2.txt
	
	
#Now convert to AFNI format
	timing_tool.py -fsl_timing_files RT_correctWrong_*.txt -write_timing RT_correctWrong_all_for_IM.1D

	cd ../..

done
