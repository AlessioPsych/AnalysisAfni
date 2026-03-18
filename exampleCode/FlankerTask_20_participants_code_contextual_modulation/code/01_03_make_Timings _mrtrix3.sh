maindir="/home/fracasso/Data/openNeuro/ds001751/derivatives/mrtrix3"

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
	cat ${subj}_task-flanker_run-1_events.tsv | awk '{if ($5=="I") {print $1, 1, 1}}' > incongruent_run1.txt
	cat ${subj}_task-flanker_run-1_events.tsv | awk '{if ($5=="C") {print $1, 1, 1}}' > congruent_run1.txt
	cat ${subj}_task-flanker_run-1_events.tsv | awk '{if ($4=="MI") {print $1, 1, 1}}' > MI_run1.txt
	cat ${subj}_task-flanker_run-1_events.tsv | awk '{if ($4=="MC") {print $1, 1, 1}}' > MC_run1.txt
	cat ${subj}_task-flanker_run-1_events.tsv | awk '{if ($5=="cue") {print $1, 1, 1}}' > cue_run1.txt
	cat ${subj}_task-flanker_run-1_events.tsv | awk '{if ($5=="blank") {print $1, 1, 1}}' > blank_run1.txt
	

	cat ${subj}_task-flanker_run-2_events.tsv | awk '{if ($5=="I") {print $1, 1, 1}}' > incongruent_run2.txt
	cat ${subj}_task-flanker_run-2_events.tsv | awk '{if ($5=="C") {print $1, 1, 1}}' > congruent_run2.txt
	cat ${subj}_task-flanker_run-2_events.tsv | awk '{if ($4=="MI") {print $1, 1, 1}}' > MI_run2.txt
	cat ${subj}_task-flanker_run-2_events.tsv | awk '{if ($4=="MC") {print $1, 1, 1}}' > MC_run2.txt
	cat ${subj}_task-flanker_run-2_events.tsv | awk '{if ($5=="cue") {print $1, 1, 1}}' > cue_run2.txt
	cat ${subj}_task-flanker_run-2_events.tsv | awk '{if ($5=="blank") {print $1, 1, 1}}' > blank_run2.txt


#Now convert to AFNI format
	timing_tool.py -fsl_timing_files congruent*.txt -write_timing congruent.1D
	timing_tool.py -fsl_timing_files incongruent*.txt -write_timing incongruent.1D
	timing_tool.py -fsl_timing_files MI_*.txt -write_timing MI.1D
	timing_tool.py -fsl_timing_files MC_*.txt -write_timing MC.1D
	timing_tool.py -fsl_timing_files cue_*.txt -write_timing cue.1D
	timing_tool.py -fsl_timing_files blank_*.txt -write_timing blank.1D

	cd ../..

done
