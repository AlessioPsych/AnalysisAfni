# use this to import DICOMs
subID <- "20230530_SDT04"
#subID <- "20230601_DLT05"

instr <- sprintf('dcm2niix_afni -o "/analyse/Project0370/%s/" /raw/Project0370/%s', subID, subID)
system( instr )

mainDir <- sprintf('/analyse/Project0370/%s', subID)
sumaDir <- sprintf('/analyse/Project0370/%s/ANATOMY/FreeSeg_result/SUMA', subID)
betaEncodingDir <- sprintf('/analyse/Project0370/%s/beta_encoding_denoised_withEndTaskPredictor.results', subID)
betaRetrievalDir <- sprintf('/analyse/Project0370/%s/beta_retrieval_denoised.results', subID)
setwd( mainDir )
getwd()

# Organise organise files in folders to have the following structure:
# 1. copy all E and R Amplitude NIFTIs into a folder called /EPI
# 2. copy all E and R Phase NIFTIs into a folder called /EPI_PHASE
# 3. copy all E and R Amp Jsons into a folder called /EPI_jsons
# 4. copy all E and R Phase Jsons into a folder called /EPI_PHASE_jsons
# 5. copy all TOPUP Amp NIFTIs into a folder called /TOPUP
# 6. copy all TOPUP Phase NIFTIs into a folder called /TOPUP_PHASE

# denoising with mrtrix
instr <- 'denoiseData.sh EPI/ *.nii 4 EPI_denoised/'; system( instr )

# slice time correction, denoised data
instr <- 'timeSliceCorrection.sh EPI_denoised/ *.nii EPI_jsons/ EPI_denoised_timeSliceCorrected/'; system( instr )

# slice time correction, original data
instr <- 'timeSliceCorrection.sh EPI/ *.nii EPI_jsons/ EPI_original_timeSliceCorrected/'; system( instr )

# compute topUp correction field
instr <- 'motionCorrect.afni.blip.sh EPI/ TOPUP/ 2 1 2-4-6 1-2-3 7 -1'; system( instr )

# denoised data: estimates motion parameters, combine them with top-up and apply motion correctoin and distortion correction in a single step
instr <- 'motionCorrectEPI.with.topUp.sh EPI_denoised_timeSliceCorrected/ topUpDir/ 2'; system( instr )
instr <- 'mv  motionCorrectEpi/ motionCorrectEpi_denoised/'; system( instr )
instr <- 'mv  motionCorrect_topUp_Epi/ motionCorrect_topUp_Epi_denoised/'; system( instr )

allFilesDenoised_BRIK <- dir( 'motionCorrect_topUp_Epi_denoised/', pattern='.BRIK' )
allFilesDenoised_HEAD <- dir( 'motionCorrect_topUp_Epi_denoised/', pattern='.HEAD' )
encodingFilesIdx <- c(1,3,5,7)
retrievalFilesIdx <- c(2,4,6,8)
dir.create( 'motionCorrect_topUp_Epi_denoised_encoding/' )
dir.create( 'motionCorrect_topUp_Epi_denoised_retrieval/' ) 
for (i in 1:length( encodingFilesIdx ) ) {
	instr <- sprintf('cp motionCorrect_topUp_Epi_denoised/%s motionCorrect_topUp_Epi_denoised_encoding/%s', allFilesDenoised_BRIK[ encodingFilesIdx[ i ] ], allFilesDenoised_BRIK[ encodingFilesIdx[ i ] ] )
	print( instr )
	system( instr )
	instr <- sprintf('cp motionCorrect_topUp_Epi_denoised/%s motionCorrect_topUp_Epi_denoised_encoding/%s', allFilesDenoised_HEAD[ encodingFilesIdx[ i ] ], allFilesDenoised_HEAD[ encodingFilesIdx[ i ] ] )
	print( instr )
	system( instr )
}

for (i in 1:length( retrievalFilesIdx ) ) {
	instr <- sprintf('cp motionCorrect_topUp_Epi_denoised/%s motionCorrect_topUp_Epi_denoised_retrieval/%s', allFilesDenoised_BRIK[ retrievalFilesIdx[ i ] ], allFilesDenoised_BRIK[ retrievalFilesIdx[ i ] ] )
	print( instr )
	system( instr )
	instr <- sprintf('cp motionCorrect_topUp_Epi_denoised/%s motionCorrect_topUp_Epi_denoised_retrieval/%s', allFilesDenoised_HEAD[ retrievalFilesIdx[ i ] ], allFilesDenoised_HEAD[ retrievalFilesIdx[ i ] ] )
	print( instr )
	system( instr )
}



# single beta estimation on denoised - for Encoding including regressors for distractor tasks ... code from Alessio ... let's see whether it runs ... seems to run!
rString01 <- paste( replicate(15, "\u0027BLOCK(5,1)\u0027"), collapse = " " )rString02 <- "\u0027BLOCK(60,1)\u0027"rString03 <- paste( replicate(15, "\u0027BLOCK(5,1)\u0027"), collapse = " " )rString04 <- "\u0027BLOCK(60,1)\u0027"
rString05 <- paste( replicate(15, "\u0027BLOCK(5,1)\u0027"), collapse = " " )rString06 <- "\u0027BLOCK(60,1)\u0027"rString07 <- paste( replicate(15, "\u0027BLOCK(5,1)\u0027"), collapse = " " )rString08 <- "\u0027BLOCK(60,1)\u0027"

instr <- paste('afni_proc.py',               '-subj_id beta_encoding_denoised_withEndTaskPredictor',               '-dsets motionCorrect_topUp_Epi_denoised_encoding/*.BRIK',               '-blocks mask scale regress',               '-regress_stim_times Stims_encoding/*.txt',               '-regress_local_times',               '-regress_polort 3',               '-regress_stim_labels T01 T02 T03 T04 T05 T06 T07 T08 T09 T10 T11 T12 T13 T14 T15 T1Dist T016 T017 T018 T019 T020 T021 T022 T023 T024 T25 T26 T27 T28 T29 T30 T2Dist T031 T032 T033 T034 T035 T036 T037 T038 T039 T40 T41 T42 T43 T44 T45 T3Dist T046 T047 T048 T049 T050 T051 T052 T053 T054 T55 T56 T57 T58 T59 T60 T4Dist',               sprintf('-regress_basis_multi %s %s %s %s %s %s %s %s', rString01, rString02, rString03, rString04, rString05, rString06, rString07, rString08),               '-regress_opts_3dD -jobs 10',               '-regress_apply_mask',               '-regress_compute_fitts',               '-regress_reml_exec',               '-regress_run_clustsim no',               '-execute')system( instr )


# single beta estimation on denoised - retrieval
instr <- paste('afni_proc.py',
		'-subj_id beta_retrieval_denoised',
		'-dsets motionCorrect_topUp_Epi_denoised_retrieval/*.BRIK',
		'-blocks mask scale regress',
		'-regress_local_times',               '-regress_polort 3',		
		'-regress_stim_times Stims_retrieval/*.txt',
		'-regress_stim_labels T01 T02 T03 T04 T05 T06 T07 T08 T09 T10 T11 T12 T13 T14 T15 T16 T17 T18 T19 T20 T21 T22 T23 T24 T25 T26 T27 T28 T29 T30 T031 T032 T033 T034 T035 T036 T037 T038 T039 T40 T41 T42 T43 T44 T45 T46 T47 T48 T49 T50 T51 T52 T53 T54 T55 T56 T57 T58 T59 T60 T061 T062 T063 T064 T065 T066 T067 T068 T069 T70 T71 T72 T73 T74 T75 T76 T77 T78 T79 T80 T81 T82 T83 T84 T85 T86 T87 T88 T89 T90 T091 T092 T093 T094 T095 T096 T097 T098 T099 T100 T101 T102 T103 T104 T105 T106 T107 T108 T109 T110 T111 T112 T113 T114 T115 T116 T117 T118 T119 T120',
		sprintf('-regress_basis \u0027BLOCK(5,1)\u0027'),
		'-regress_opts_3dD -jobs 10',
		'-regress_apply_mask',
		'-regress_compute_fitts',
		'-regress_reml_exec',
		'-regress_run_clustsim no',
		'-execute')
system( instr )


# run this
# one beta for visual stimulus onset encoding
instr <- paste('afni_proc.py',
		'-subj_id beta_vis_enc_denoised',
		'-dsets motionCorrect_topUp_Epi_denoised_encoding/*.BRIK',
		'-blocks mask scale regress',
		'-regress_stim_times Stims_vis_encoding/*.txt',
		'-regress_stim_labels P_onset',
		sprintf('-regress_basis \u0027BLOCK(5,1)\u0027'),
		'-regress_opts_3dD -jobs 10',
		'-regress_apply_mask',
		'-regress_compute_fitts',
		'-regress_reml_exec',
		'-regress_run_clustsim no',
		'-execute')
system( instr )

# one beta for visual stimulus onset retrieval
instr <- paste('afni_proc.py',
		'-subj_id beta_vis_retr_denoised',
		'-dsets motionCorrect_topUp_Epi_denoised_retrieval/*.BRIK',
		'-blocks mask scale regress',
		'-regress_stim_times Stims_vis_retrieval/*.txt',
		'-regress_stim_labels P_onset',
		sprintf('-regress_basis \u0027BLOCK(5,1)\u0027'),
		'-regress_opts_3dD -jobs 10',
		'-regress_apply_mask',
		'-regress_compute_fitts',
		'-regress_reml_exec',
		'-regress_run_clustsim no',
		'-execute')
system( instr )


# Set up a simple memory contrast between remembered and forgotten (collapsed accross encoding and retrieval) 
instr <- paste('afni_proc.py',
		'-subj_id beta_memory',
		'-dsets motionCorrect_topUp_Epi_original/*.BRIK',
		'-blocks mask scale regress',
		'-regress_stim_times Stims_Memory/*.txt',
		'-regress_stim_labels rem_pred forg_pred',
		sprintf('-regress_basis \u0027BLOCK(5,1)\u0027'),
		'-regress_opts_3dD -gltsym \u0027SYM: +rem_pred -forg_pred\u0027 -glt_label 1 Hits_vs_Miss -jobs 10',              
                '-regress_apply_mask',
		'-regress_compute_fitts',
		'-regress_reml_exec',
		'-regress_run_clustsim no',
		'-execute')
system( instr )


# coregistration: see folder Coregistration/ and script coregistrationScript.sh

# bring betas (encoding) into anatomy (suma folder) ... # Changed output directory to base directory because I don't have write access to SUMA Dir
instr <- sprintf( '3dAllineate -prefix %s/statsEncoding.nii.gz -1Dmatrix_apply Coregistration/coregMat.1D -final NN -input %s/stats.beta_encoding_denoised_withEndTaskPredictor+orig -master %s/anatCopy.nii.gz', sumaDir, betaEncodingDir, sumaDir )
system( instr )

# bring betas (retrieval) into anatomy (suma folder) ...  # Changed output directory to base directory because I don't have write access to SUMA Dir
instr <- sprintf( '3dAllineate -prefix %s/statsRetrieval.nii.gz -1Dmatrix_apply Coregistration/coregMat.1D -final NN -input %s/stats.beta_retrieval_denoised+orig -master %s/anatCopy.nii.gz', sumaDir, betaRetrievalDir, sumaDir )
system( instr )

# -- from here on it's Matlab







#####################################################
########### OLD STUFF THAT IS NOT NEEDed ############
#####################################################




# Coregistration
# first average EPIs to get a better signal to noise ratio - we do this here for all preprocessed EPIs - Encoding and Retrieval together
# this needs to be run in Bash so make sure you're outta R
computeAmplitudeAnatomy.sh /analyse/Project0370/20230418_CME20/motionCorrect_topUp_Epi_denoised/ amplitudeAnatomy.nii
# Copy EPIamplitudes and segmented brain.nii in same folder 
mkdir /analyse/Project0370/20230418_CME20/Coregistration
cp /analyse/Project0370/20230418_CME20/amplitudeAnatomy.nii /analyse/Project0370/20230418_CME20/Coregistration/
# This step is manual open AFNI and use -> Define Datamode -> PLugin -> Nudge to manually aling the EPIs to the Brain and write resulting file in Co-registration - you'll get a warning but that's fine, we sill have the original file in the other folder
instr <- paste('3dAllineate',
	'-base brain.nii',
	'-source amplitudeAnatomy.nii',
	'-onepass',
	'-warp sho',	
	'-prefix amplitudeAnatomy_cor2.nii',
	'-1Dmatrix_save matrix.1D',
	'-cmass',
	'-maxrot 4',
	' -maxshf 4',
	'float')
system( instr )

	  

5


# Retired stuff and things that didn't work down here ....




# single beta estimation try denoised - solve collinearity problem by running it separately for Encoding ... 
instr <- paste('afni_proc.py', # basic AFNI command to run the analysis
		'-subj_id beta_all_encoding_denoised', # specify output directory where results are stored
		'-dsets motionCorrect_topUp_Epi_denoised_encoding/*.BRIK', # specify input data files; the wildcard indicates multiple runs are provided as input; AFNI will read these in in alphabetical order; make sure your regressors match these, and make sure to appropriately set regress_local_times
		'-blocks mask scale regress', # this specifies the processing steps to be performaed; The 'mask' block creates a brain mask to exclude non-brain voxles, the 'scale' block scales the time series to a baseline, and the 'regress' block performs regression analysis;  
		'-regress_stim_times Stims_encoding/*.txt', # stimulus timings stored in txt files; 
		'-regress_local_times', # tells AFNI to use local times, ie. interpret timings in *.txt files with respect to start of run;
		'-regress_polort 3', # perform polynomial fit (3 = cubic) to get rid of low freq fluctuations
		'-regress_stim_labels T1 T2 T3 T4 T5 T6 T7 T8 T9 T10 T11 T12 T13 T14 T15 T16 T17 T18 T19 T20 T21 T22 T23 T24 T25 T26 T27 T28 T29 T30', # stimulus labels
		sprintf('-regress_basis \u0027BLOCK(5,1)\u0027'), # this sets the basis function for the regression analysis; The BLOCK(5,1) specifies a block design with a duration of 5 seconds and a delay of  second
		'-regress_opts_3dD -jobs 10', # Enables additional regression analysis options and job -10 specifies the number of parallel processing jobs to use
		'-regress_apply_mask', # applies the brain mask created in the 'mask' block to restrict analysis to brain vxls
		'-regress_compute_fitts', # comput goodness of fit stats
		'-regress_reml_exec', # performs a restricteds max likleihood (REML) estimatation to obtain model parameter estimates; do we need that?
		'-regress_run_clustsim no', # disables cluster size simulation, which is used to estimate the spatial extent of activation clusters
		'-execute') # executes the processing steps above
system( instr )



# run this
# single beta estimation try, original
instr <- paste('afni_proc.py',
		'-subj_id betaAll_original',
		'-dsets motionCorrect_topUp_Epi_original/*.BRIK',
		'-blocks mask scale regress',
		'-regress_stim_times Stims/*.txt',
		'-regress_stim_labels T01 T02 T03 T04 T05',
		sprintf('-regress_basis \u0027BLOCK(5,1)\u0027'),
		'-regress_opts_3dD -jobs 10',
		'-regress_apply_mask',
		'-regress_compute_fitts',
		'-regress_reml_exec',
		'-regress_run_clustsim no',
		'-execute')
system( instr )

# run this
# one beta for visual stimulus onset (encoding and retrieval)
instr <- paste('afni_proc.py',
		'-subj_id beta_vis_denoised',
		'-dsets motionCorrect_topUp_Epi_denoised/*.BRIK',
		'-blocks mask scale regress',
		'-regress_stim_times Stims_all4one/*.txt',
		'-regress_stim_labels P_onset',
		sprintf('-regress_basis \u0027BLOCK(5,1)\u0027'),
		'-regress_opts_3dD -jobs 10',
		'-regress_apply_mask',
		'-regress_compute_fitts',
		'-regress_reml_exec',
		'-regress_run_clustsim no',
		'-execute')
system( instr )

# run this
# one beta for visual stimulus onset on original data (encoding and retrieval)
instr <- paste('afni_proc.py',
		'-subj_id beta_vis_original',
		'-dsets motionCorrect_topUp_Epi_original/*.BRIK',
		'-blocks mask scale regress',
		'-regress_stim_times Stims_all4one/*.txt',
		'-regress_stim_labels P_onset',
		sprintf('-regress_basis \u0027BLOCK(5,1)\u0027'),
		'-regress_opts_3dD -jobs 10',
		'-regress_apply_mask',
		'-regress_compute_fitts',
		'-regress_reml_exec',
		'-regress_run_clustsim no',
		'-execute')
system( instr )



