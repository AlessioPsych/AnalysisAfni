rm( list=ls() ); gc();

#### Move to the desired folder ####

sprintf('Move to the desired folder:')
homeDir <- '/scratch/af4887/Proj_Gaia_David'
mainDir <- '/scratch/af4887/Proj_Gaia_David/afni_processed'
inputFreesurfer <- '/scratch/af4887/Proj_Gaia_David/freesurfer'

setwd( homeDir )
sprintf( 'Current folder: %s', print( getwd() ) )
runCodeFlag <- 1

if ( dir.exists('afni_glm') ) {
  instr <- 'rm -R afni_glm/'
  print( instr )
  if ( runCodeFlag==1 ) { system( instr ) }
}
instr <- 'mkdir afni_glm/'
print( instr )
if ( runCodeFlag==1 ) { system( instr ) }

setwd( mainDir )
singleSubjectFolders <- list.dirs( recursive=FALSE )

setwd( inputFreesurfer )
singleSubjectFolders_freesurfer <- list.dirs( recursive=FALSE )

args <- commandArgs(trailingOnly = TRUE)
print(args[1])
var<- args[1]

for ( nSubj in var : var  ) { # nSubj <- 1 length( singleSubjectFolders )
  
  nSubj <- 1
  
  # convert brain.mgz from freesurfer folder
  setwd( sprintf( '%s/%s/mri', inputFreesurfer, singleSubjectFolders_freesurfer[ nSubj ] ) ) 
  print( getwd() )
  if ( file.exists( 'brain.nii' ) ) {
    print( sprintf('file brain.nii exists, remove...') )
    if (runCodeFlag==1) { system( 'rm brain.nii' ) }
  }
  instr <- 'mri_convert brain.mgz brain.nii'
  print( instr )
  if (runCodeFlag==1) { system( instr ) }

  setwd( homeDir )
  setwd( 'afni_glm' )
  splitName <- strsplit( singleSubjectFolders[ nSubj ], '[./]' )[[1]][3]
  subjName <- sprintf('%s-glm', splitName )
  
  instr <- paste('afni afni_proc.py',
                 sprintf( '-subj_id %s', subjName ),
                 #sprintf( '-dsets /scratch/af4887/Proj_Gaia_David/afni_processed/sub-03.results/pb03.sub-03.r01.volreg+orig.BRIK /scratch/af4887/Proj_Gaia_David/afni_processed/sub-03.results/pb03.sub-03.r02.volreg+orig.BRIK' ),
                 sprintf( '-dsets %s/%s/*.volreg+orig.BRIK', mainDir, singleSubjectFolders[ nSubj ] ),
                 sprintf('-copy_anat %s/%s/mri/brain.nii', inputFreesurfer, singleSubjectFolders_freesurfer[ nSubj ] ),
                 #'-blocks align blur mask scale regress',
                 '-blocks blur scale regress',
                 #'-align_opts_aea -epi2anat',
                 '-blur_size 2.5',
                 sprintf('-regress_stim_times %s/%s/*.1D', homeDir, 'del_stim_times'),
                 '-regress_local_times',
                 '-regress_polort 5',
                 '-regress_stim_types AM2',
                 '-regress_stim_labels T01 T02',
                 sprintf('-regress_basis \u0027dmUBLOCK\u0027'),
                 '-regress_opts_3dD -jobs 1',
                 #'-mask_apply epi',
                 #'-regress_compute_fitts',
                 '-regress_run_clustsim no',
                 '-execute')
  print( instr )
  system( instr )
  
  # rm brain.nii from freesurfer folder
  setwd( sprintf( '%s/%s/mri', inputFreesurfer, singleSubjectFolders_freesurfer[ nSubj ] ) ) 
  print( getwd() )
  if ( file.exists( 'brain.nii' ) ) {
    print( sprintf('file brain.nii exists, remove...') )
    if (runCodeFlag==1) { system( 'rm brain.nii' ) }
  }
  
} 


