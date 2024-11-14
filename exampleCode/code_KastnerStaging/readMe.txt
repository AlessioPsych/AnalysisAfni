01_00_processingEPI_afni.R

	runs denoise with mritrix, timeslice correction, distortion correction and motion correction on each session. Saves corrected data and motion parameters
	
01_01_processingEPI_afni_extraData.R	

	runs denoise with mritrix, timeslice correction, distortion correction and motion correction on each session. Saves corrected data and motion parameters. On extra data (eye 		bar 2 and eyebar 3).

02_00_coregistration.R

	estimates the coregitration matrix of epis on individual anatomy;

02_01_coregistration_extraData.R

	estimates the coregitration matrix of epis on individual anatomy; On extra data (eyebar2 and eyebar3)	

03_00_bringEPIOnAnat.R

	interpolates (ilnear interpolation), each epi on the individual anatomy at a resolution of 1.5mm isotropic (the original resolution)

03_01_bringEPIOnAnat_extraData.R	

	interpolates (ilnear interpolation), each epi on the individual anatomy at a resolution of 1.5mm isotropic (the original resolution), on extra data (eyebar2 and  eyebar3)	

04_00_averageExpConditions.R (ONLY FOR LIMITED CONDITIONS, DO NOT USE)

	!!! here concatenate missing data from a couple of runs, specifically on session 21 and session 26. See file !!!!

	averages all individuals that took part in each experimental condition in the dataset and stores the average in a newly created folder under 'expConditions'

04_01_averageExpConditions.R (FOR ALL CONDITIONS, USE THIS ONE)

	!!! here concatenate missing data from a couple of runs, specifically on session 21 and session 26. See file !!!!

	averages all individuals that took part in each experimental condition in the dataset and stores the average in a newly created folder under 'expConditions'

05_KastnerClassicAnalysis.R 

	performs kastner classic analysis and saves volumetric data
	# to do: invert phase fit outcome to match the modelling results!

06_KastnerClassicBringVolumeOnSurface.R

	copy volumetric data output into individual subject SUMA folder and convert into surfaces (niml and 1D)

07_KastnerClassicAverageOnSurface.R

	computes averages on standard surface

08_KastnerClassicGetRoisOnIndividualParticipants.R

	computes volumetric ROIS from Wang2015 atlas on individual participants anatomy folder
	
09_KastnerClassicGetParticipantsROIData

	stores individual voxel modelling results, per ROI (from Wang2015), from all Kastner Classic models, including Cw and 		CCw, on separate RData files
	
10_KastnerClassic_pRF_Analysis.R

	runs pRF analysis on all participants where pRF data is available
	
11_pRFBringVolumeOnSurface.R

	copy pRF volumetric data output into individual subject SUMA folder and convert into surfaces (niml and 1D)

12_pRFAverageOnSurface.R

	computes pRF averages on standard surface
	
13_pRFGetRoisOnIndividualParticipants.R

	computes volumetric ROIS from Wang2015 atlas on individual participants anatomy folder

14_pRFGetParticipantsROIData.R

	stores individual voxel modelling results, per ROI (from Wang2015), from pRF models, on a single RDATA file in mainFolder/results
	
15_analysis_rois_Kastner_classic.R

	perform analysis for Kastner Classic paper	

