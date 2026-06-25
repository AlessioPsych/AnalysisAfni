rm(list=ls())

subID <- "20230620_CCY20"

# set directories
mainDir <- sprintf('/analyse/Project0370/%s', subID)
sumaDir <- sprintf('/analyse/Project0370/%s/ANATOMY/FreeSeg_result/SUMA', subID)
betaEncodingDir <- sprintf('/analyse/Project0370/%s/beta_encoding_phase_denoised_withEndTaskPredictor.results', subID)
betaRetrievalDir <- sprintf('/analyse/Project0370/%s/beta_retrieval_phase_denoised.results', subID)
setwd( mainDir )
getwd()

# clean up 
system('rm proc.beta_encoding_phase_denoised_withEndTaskPredictor')
system('rm output.proc.beta_encoding_phase_denoised_withEndTaskPredictor')
system('rm proc.beta_retrieval_phase_denoised')
system('rm output.proc.beta_retrieval_phase_denoised')
system('rm proc.beta_vis_enc_phase_denoised')
system('rm output.proc.beta_vis_enc_phase_denoised')
system('rm proc.beta_vis_retr_phase_denoised')
system('rm output.proc.beta_vis_retr_phase_denoised')
system('rm -R EPI_phase_denoised/')
system('rm -R EPI_phase_denoised_timeSliceCorrected')
system('rm -R motionCorrectEpi_phase_denoised/')
system('rm -R motionCorrect_topUp_Epi_phase_denoised')
system('rm -R motionCorrect_topUp_Epi_phase_denoised_encoding')
system('rm -R motionCorrect_topUp_Epi_phase_denoised_retrieval')
system('rm -R beta_encoding_phase_denoised_withEndTaskPredictor.results')
system('rm -R beta_retrieval_phase_denoised.results')
system('rm -R beta_vis_enc_phase_denoised.results')
system('rm -R beta_vis_retr_phase_denoised.results')
system( sprintf( 'rm %s/statsEncoding_phaseDenoised.nii.gz', sumaDir ) )
system( sprintf( 'rm %s/statsRetrieval_phaseDenoised.nii.gz', sumaDir ) )
system( sprintf( 'rm %s/statsEncoding_phaseDenoised_REML.nii.gz', sumaDir ) )
system( sprintf( 'rm %s/statsRetrieval_phaseDenoised_REML.nii.gz', sumaDir ) )

# denoising with mrtrix
instr <- 'denoiseData_with_phase.sh EPI/ EPI_PHASE/ *.nii 8 EPI_phase_denoised/'
print( instr )
system( instr )

# slice time correction, denoised data
instr <- 'timeSliceCorrection_noRJSON.sh EPI_phase_denoised/ *.nii EPI_jsons/ EPI_phase_denoised_timeSliceCorrected/' 
print( instr )
system( instr )

# compute topUp correction field (already done with the previous analysis, no need to do it again)
# we already have topUpDir from the previous analysis, we can use it here.
#instr <- 'motionCorrect.afni.blip.sh EPI/ TOPUP/ 2 1 2-4-6 1-2-3 7 -1'; system( instr )

# denoised data: estimates motion parameters, combine them with top-up and apply motion correctoin and distortion correction in a single step
instr <- 'motionCorrectEPI.with.topUp.sh EPI_phase_denoised_timeSliceCorrected/ topUpDir/ 2' 
print( instr )
system( instr )
instr <- 'mv  motionCorrectEpi/ motionCorrectEpi_phase_denoised/'
print( instr )
system( instr )
instr <- 'mv  motionCorrect_topUp_Epi/ motionCorrect_topUp_Epi_phase_denoised/'
print( instr )
system( instr )

allFilesDenoised_BRIK <- dir( 'motionCorrect_topUp_Epi_phase_denoised/', pattern='.BRIK' )
allFilesDenoised_HEAD <- dir( 'motionCorrect_topUp_Epi_phase_denoised/', pattern='.HEAD' )
encodingFilesIdx <- c(1,3,5,7)
retrievalFilesIdx <- c(2,4,6,8)
dir.create( 'motionCorrect_topUp_Epi_phase_denoised_encoding/' )
dir.create( 'motionCorrect_topUp_Epi_phase_denoised_retrieval/' ) 
for (i in 1:length( encodingFilesIdx ) ) {
	instr <- sprintf('cp motionCorrect_topUp_Epi_phase_denoised/%s motionCorrect_topUp_Epi_phase_denoised_encoding/%s', allFilesDenoised_BRIK[ encodingFilesIdx[ i ] ], allFilesDenoised_BRIK[ encodingFilesIdx[ i ] ] )
	print( instr )
	system( instr )
	instr <- sprintf('cp motionCorrect_topUp_Epi_phase_denoised/%s motionCorrect_topUp_Epi_phase_denoised_encoding/%s', allFilesDenoised_HEAD[ encodingFilesIdx[ i ] ], allFilesDenoised_HEAD[ encodingFilesIdx[ i ] ] )
	print( instr )
	system( instr )
}

for (i in 1:length( retrievalFilesIdx ) ) {
	instr <- sprintf('cp motionCorrect_topUp_Epi_phase_denoised/%s motionCorrect_topUp_Epi_phase_denoised_retrieval/%s', allFilesDenoised_BRIK[ retrievalFilesIdx[ i ] ], allFilesDenoised_BRIK[ retrievalFilesIdx[ i ] ] )
	print( instr )
	system( instr )
	instr <- sprintf('cp motionCorrect_topUp_Epi_phase_denoised/%s motionCorrect_topUp_Epi_phase_denoised_retrieval/%s', allFilesDenoised_HEAD[ retrievalFilesIdx[ i ] ], allFilesDenoised_HEAD[ retrievalFilesIdx[ i ] ] )
	print( instr )
	system( instr )
}

# single beta estimation on phase denoised - for Encoding including regressors for distractor tasks ... code from Alessio ... let's see whether it runs ... seems to run!
rString01 <- paste( replicate(15, "\u0027BLOCK(5,1)\u0027"), collapse = " " )
rString02 <- "\u0027BLOCK(60,1)\u0027"
rString03 <- paste( replicate(15, "\u0027BLOCK(5,1)\u0027"), collapse = " " )
rString04 <- "\u0027BLOCK(60,1)\u0027"
rString05 <- paste( replicate(15, "\u0027BLOCK(5,1)\u0027"), collapse = " " )
rString06 <- "\u0027BLOCK(60,1)\u0027"
rString07 <- paste( replicate(15, "\u0027BLOCK(5,1)\u0027"), collapse = " " )
rString08 <- "\u0027BLOCK(60,1)\u0027"

instr <- paste('afni_proc.py',
               '-subj_id beta_encoding_phase_denoised_withEndTaskPredictor',
               '-dsets motionCorrect_topUp_Epi_phase_denoised_encoding/*.BRIK',
               '-blocks mask scale regress',
               '-regress_stim_times Stims_encoding/*.txt',
               '-regress_local_times',
               '-regress_polort 3',
               '-regress_stim_labels T01 T02 T03 T04 T05 T06 T07 T08 T09 T10 T11 T12 T13 T14 T15 T1Dist T016 T017 T018 T019 T020 T021 T022 T023 T024 T25 T26 T27 T28 T29 T30 T2Dist T031 T032 T033 T034 T035 T036 T037 T038 T039 T40 T41 T42 T43 T44 T45 T3Dist T046 T047 T048 T049 T050 T051 T052 T053 T054 T55 T56 T57 T58 T59 T60 T4Dist',
               sprintf('-regress_basis_multi %s %s %s %s %s %s %s %s', rString01, rString02, rString03, rString04, rString05, rString06, rString07, rString08),
               '-regress_opts_3dD -jobs 10',
               '-regress_apply_mask',
               '-regress_compute_fitts',
               '-regress_reml_exec',
               '-regress_run_clustsim no',
               '-execute')
system( instr )

# single beta estimation on denoised - retrieval
instr <- paste('afni_proc.py',
		'-subj_id beta_retrieval_phase_denoised',
		'-dsets motionCorrect_topUp_Epi_phase_denoised_retrieval/*.BRIK',
		'-blocks mask scale regress',
		'-regress_local_times',
               '-regress_polort 3',		
		'-regress_stim_times Stims_retrieval/*.txt',
		'-regress_stim_labels T01 T02 T03 T04 T05 T06 T07 T08 T09 T10 T11 T12 T13 T14 T15 T16 T17 T18 T19 T20 T21 T22 T23 T24 T25 T26 T27 T28 T29 T30 T031 T032 T033 T034 T035 T036 T037 T038 T039 T40 T41 T42 T43 T44 T45 T46 T47 T48 T49 T50 T51 T52 T53 T54 T55 T56 T57 T58 T59 T60 T061 T062 T063 T064 T065 T066 T067 T068 T069 T70 T71 T72 T73 T74 T75 T76 T77 T78 T79 T80 T81 T82 T83 T84 T85 T86 T87 T88 T89 T90 T091 T092 T093 T094 T095 T096 T097 T098 T099 T100 T101 T102 T103 T104 T105 T106 T107 T108 T109 T110 T111 T112 T113 T114 T115 T116 T117 T118 T119 T120',
		sprintf('-regress_basis \u0027BLOCK(5,1)\u0027'),
		'-regress_opts_3dD -jobs 10',
		'-regress_apply_mask',
		'-regress_compute_fitts',
		'-regress_reml_exec',
		'-regress_run_clustsim no',
		'-execute')
print( instr )
system( instr )

# run this
# one beta for visual stimulus onset encoding
instr <- paste('afni_proc.py',
		'-subj_id beta_vis_enc_phase_denoised',
		'-dsets motionCorrect_topUp_Epi_phase_denoised_encoding/*.BRIK',
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
print( instr )
system( instr )

# one beta for visual stimulus onset retrieval
instr <- paste('afni_proc.py',
		'-subj_id beta_vis_retr_phase_denoised',
		'-dsets motionCorrect_topUp_Epi_phase_denoised_retrieval/*.BRIK',
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
print( instr )
system( instr )


# Set up a simple memory contrast between remembered and forgotten (collapsed accross encoding and retrieval) 
#instr <- paste('afni_proc.py',
#		'-subj_id beta_memory',
#		'-dsets motionCorrect_topUp_Epi_original/*.BRIK',
#		'-blocks mask scale regress',
#		'-regress_stim_times Stims_Memory/*.txt',
#		'-regress_stim_labels rem_pred forg_pred',
#		sprintf('-regress_basis \u0027BLOCK(5,1)\u0027'),
#		'-regress_opts_3dD -gltsym \u0027SYM: +rem_pred -forg_pred\u0027 -glt_label 1 Hits_vs_Miss -jobs 10',              
#                '-regress_apply_mask',
#		'-regress_compute_fitts',
#		'-regress_reml_exec',
#		'-regress_run_clustsim no',
#		'-execute')
#system( instr )


# coregistration: see folder Coregistration/ and script coregistrationScript.sh

# bring betas (encoding) into anatomy (suma folder) ... # Changed output directory to base directory because I don't have write access to SUMA Dir
instr <- sprintf( '3dAllineate -prefix %s/statsEncoding_phaseDenoised.nii.gz -1Dmatrix_apply Coregistration/coregMat.1D -final NN -input %s/stats.beta_encoding_phase_denoised_withEndTaskPredictor+orig -master %s/anatCopy.nii.gz', sumaDir, betaEncodingDir, sumaDir )
system( instr )

# bring betas (retrieval) into anatomy (suma folder) ...  # Changed output directory to base directory because I don't have write access to SUMA Dir
instr <- sprintf( '3dAllineate -prefix %s/statsRetrieval_phaseDenoised.nii.gz -1Dmatrix_apply Coregistration/coregMat.1D -final NN -input %s/stats.beta_retrieval_phase_denoised+orig -master %s/anatCopy.nii.gz', sumaDir, betaRetrievalDir, sumaDir )
system( instr )

# bring betas (encoding REML) into anatomy (suma folder) ... # Changed output directory to base directory because I don't have write access to SUMA Dir
instr <- sprintf( '3dAllineate -prefix %s/statsEncoding_phaseDenoised_REML.nii.gz -1Dmatrix_apply Coregistration/coregMat.1D -final NN -input %s/stats.beta_encoding_phase_denoised_withEndTaskPredictor_REML+orig -master %s/anatCopy.nii.gz', sumaDir, betaEncodingDir, sumaDir )
system( instr )

# bring betas (retrieval REML) into anatomy (suma folder) ...  # Changed output directory to base directory because I don't have write access to SUMA Dir
instr <- sprintf( '3dAllineate -prefix %s/statsRetrieval_phaseDenoised_REML.nii.gz -1Dmatrix_apply Coregistration/coregMat.1D -final NN -input %s/stats.beta_retrieval_phase_denoised_REML+orig -master %s/anatCopy.nii.gz', sumaDir, betaRetrievalDir, sumaDir )
system( instr )

# open up for modifications:
system('chmod 777 -R EPI_phase_denoised/')
system('chmod 777 -R EPI_phase_denoised_timeSliceCorrected')
system('chmod 777 -R motionCorrectEpi_phase_denoised/')
system('chmod 777 -R motionCorrect_topUp_Epi_phase_denoised')
system('chmod 777 -R motionCorrect_topUp_Epi_phase_denoised_encoding')
system('chmod 777 -R motionCorrect_topUp_Epi_phase_denoised_retrieval')
system('chmod 777 -R beta_encoding_phase_denoised_withEndTaskPredictor.results')
system('chmod 777 -R beta_retrieval_phase_denoised.results')
system('chmod 777 -R beta_vis_enc_phase_denoised.results')
system('chmod 777 -R beta_vis_retr_phase_denoised.results')
system( sprintf( 'chmod 777 %s/statsEncoding_phaseDenoised.nii.gz', sumaDir ) )
system( sprintf( 'chmod 777 %s/statsRetrieval_phaseDenoised.nii.gz', sumaDir ) )
system( sprintf( 'chmod 777 %s/statsEncoding_phaseDenoised_REML.nii.gz', sumaDir ) )
system( sprintf( 'chmod 777 %s/statsRetrieval_phaseDenoised_REML.nii.gz', sumaDir ) )
system('chmod 777 proc.beta_encoding_phase_denoised_withEndTaskPredictor')
system('chmod 777 output.proc.beta_encoding_phase_denoised_withEndTaskPredictor')
system('chmod 777 proc.beta_retrieval_phase_denoised')
system('chmod 777 output.proc.beta_retrieval_phase_denoised')
system('chmod 777 proc.beta_vis_enc_phase_denoised')
system('chmod 777 output.proc.beta_vis_enc_phase_denoised')
system('chmod 777 proc.beta_vis_retr_phase_denoised')
system('chmod 777 output.proc.beta_vis_retr_phase_denoised')
# -- from here on it's Matlab
