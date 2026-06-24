maindir="/mnt/disk01/ds001751_FlankerTask_context/derivatives/mrtrix3"

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
	cat ${subj}_task-flanker_run-1_events.tsv | awk '{if ($5=="I" && $4=="MC" && $7=="1") {print $1, 1, 1}}' > incongruent_MC_run1.txt
	cat ${subj}_task-flanker_run-1_events.tsv | awk '{if ($5=="C" && $4=="MC" && $7=="1") {print $1, 1, 1}}' > congruent_MC_run1.txt
	cat ${subj}_task-flanker_run-1_events.tsv | awk '{if ($5=="I" && $4=="MI" && $7=="1") {print $1, 1, 1}}' > incongruent_MI_run1.txt
	cat ${subj}_task-flanker_run-1_events.tsv | awk '{if ($5=="C" && $4=="MI" && $7=="1") {print $1, 1, 1}}' > congruent_MI_run1.txt
	#cat ${subj}_task-flanker_run-1_events.tsv | awk '{if ($4=="MI") {print $1, 1, 1}}' > MI_run1.txt
	#cat ${subj}_task-flanker_run-1_events.tsv | awk '{if ($4=="MC") {print $1, 1, 1}}' > MC_run1.txt
	cat ${subj}_task-flanker_run-1_events.tsv | awk '{if ($5=="cue") {print $1, 1, 1}}' > cue_run1.txt
	cat ${subj}_task-flanker_run-1_events.tsv | awk '{if ($5=="blank") {print $1, 1, 1}}' > blank_run1.txt
	

	cat ${subj}_task-flanker_run-2_events.tsv | awk '{if ($5=="I" && $4=="MC" && $7=="1") {print $1, 1, 1}}' > incongruent_MC_run2.txt
	cat ${subj}_task-flanker_run-2_events.tsv | awk '{if ($5=="C" && $4=="MC" && $7=="1") {print $1, 1, 1}}' > congruent_MC_run2.txt
	cat ${subj}_task-flanker_run-2_events.tsv | awk '{if ($5=="I" && $4=="MI" && $7=="1") {print $1, 1, 1}}' > incongruent_MI_run2.txt
	cat ${subj}_task-flanker_run-2_events.tsv | awk '{if ($5=="C" && $4=="MI" && $7=="1") {print $1, 1, 1}}' > congruent_MI_run2.txt
	#cat ${subj}_task-flanker_run-2_events.tsv | awk '{if ($4=="MI") {print $1, 1, 1}}' > MI_run2.txt
	#cat ${subj}_task-flanker_run-2_events.tsv | awk '{if ($4=="MC") {print $1, 1, 1}}' > MC_run2.txt
	cat ${subj}_task-flanker_run-2_events.tsv | awk '{if ($5=="cue") {print $1, 1, 1}}' > cue_run2.txt
	cat ${subj}_task-flanker_run-2_events.tsv | awk '{if ($5=="blank") {print $1, 1, 1}}' > blank_run2.txt


#Now convert to AFNI format
	timing_tool.py -fsl_timing_files congruent_MC*.txt -write_timing congruent_MC.1D
	timing_tool.py -fsl_timing_files incongruent_MC*.txt -write_timing incongruent_MC.1D
	timing_tool.py -fsl_timing_files congruent_MI*.txt -write_timing congruent_MI.1D
	timing_tool.py -fsl_timing_files incongruent_MI*.txt -write_timing incongruent_MI.1D
	#timing_tool.py -fsl_timing_files MI_*.txt -write_timing MI.1D
	#timing_tool.py -fsl_timing_files MC_*.txt -write_timing MC.1D
	timing_tool.py -fsl_timing_files cue_*.txt -write_timing cue.1D
	timing_tool.py -fsl_timing_files blank_*.txt -write_timing blank.1D

	cd ../..

done
