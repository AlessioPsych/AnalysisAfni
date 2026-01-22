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

	#cat ${subj}_task-flanker_run-1_events.tsv | awk '{if ($5=="correct") {print $1, $2, 1}}' > RT_run1.txt

	#cat ${subj}_task-flanker_run-2_events.tsv | awk '{if ($5=="correct") {print $1, $2, 1}}' > RT_run2.txt
	
	cat ${subj}_task-flanker_run-1_events.tsv | awk '{if ($5=="correct") {print $1, $4, 1}}' > RT_run1.txt
	cat ${subj}_task-flanker_run-1_events.tsv | awk '{if ($5=="incorrect") {print $1, $4, 2}}' > incorrect_RT_run1.txt


	cat ${subj}_task-flanker_run-2_events.tsv | awk '{if ($5=="correct") {print $1, $4, 1}}' > RT_run2.txt
	cat ${subj}_task-flanker_run-2_events.tsv | awk '{if ($5=="incorrect") {print $1, $4, 2}}' > incorrect_RT_run2.txt

	
#Now convert to AFNI format
	rm RT_correct.1D
	rm incorrect_RT.1D
	timing_tool.py -fsl_timing_files RT_run*.txt -write_as_married -write_timing RT_correct.1D
	timing_tool.py -fsl_timing_files incorrect_RT_run*.txt -write_as_married -write_timing incorrect_RT.1D


	cd ../..

done
