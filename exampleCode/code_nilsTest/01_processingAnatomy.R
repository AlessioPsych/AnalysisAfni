rm( list=ls() )

subjNumber <- 1

if (subjNumber==1) {
  mainDir <- '/analyse/Project0226/tests/nilsTest/data/20230622_SHI27/ANATOMY'
  secondInv <- '9_mp2rage_sag_p3_0.6mm.nii.gz'
  anatomy <- '10_mp2rage_sag_p3_0.6mm.nii.gz'
}

setwd( mainDir )
print('################################################')
print('################################################')
print('################################################')
print( getwd() )
print('################################################')
print('################################################')
print('################################################')

# clean up resamped volume ad freesurfer segmentation
if ( file.exists( 'anatCopy.nii.gz' ) ) { system('rm anatCopy.nii.gz') }
if ( file.exists( 'fsaverage' ) ) { system('rm fsaverage') }
if ( dir.exists( 'anat_test' ) ) { system('rm -R anat_test') }

# get t1w file
participantFiles <- getwd()
participantFiles

#skul strip
print( '...' )
instr <- sprintf( 'skullStrip_mp2rage.sh %s %s %s', anatomy, secondInv, '0.7' ); print( instr ); system( instr )
print( '...' )

#redefine SUBJECTS_DIR to current directory
currentDir <- getwd()
print( '...' )
Sys.setenv( SUBJECTS_DIR = paste( currentDir ) )
print( Sys.getenv( 'SUBJECTS_DIR') )
print( '...' )

#recon-all instruction
print( '...' )
instr <- sprintf( 'recon-all -subjid anat_test -i anatCopy.nii.gz -all -parallel -openmp 8' ); print( instr ); system( instr )
print( '...' )

#suma instruction
print( '...' )
instr <- sprintf( '@SUMA_Make_Spec_FS -NIFTI -fspath %s/anat_test -sid anat_test', currentDir ); print( instr ); system( instr )
print( '...' )




