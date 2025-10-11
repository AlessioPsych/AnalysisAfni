rm( list=ls() )
mainDir <- '/analyse/Project0165/Fracasso_Anatomy/AHEAD_database/10007840'
partDir <- 'part1/'
volumetricCodeCall <- 'sh /analyse/Project0165/Fracasso_Anatomy/AHEAD_database/10007840/analysis_15032022/volumetric_code.sh'
profileCodeCall <- 'sh /analyse/Project0165/Fracasso_Anatomy/AHEAD_database/10007840/analysis_15032022/profile_code.sh'

setwd( mainDir )
setwd( partDir )
getwd()
participantsFolders <- dir()
for ( participantsFoldersIndex in 2:4 ) { #seq( 1,length( participantsFolders ) thickness and curvature
  
  setwd( mainDir )
  setwd( partDir )
  print( '...' )
  print( '...' )
  print( sprintf( '%s...', participantsFolders[ participantsFoldersIndex ] ) )
  print( '...' )
  print( '...' )
  setwd( sprintf( '%s', participantsFolders[ participantsFoldersIndex ] ) )
  setwd( dir()[1] )
  setwd( dir()[1] )
  print( dir() )
  setwd( 'AHEAD_test' )
  setwd( 'SUMA' )
  print( dir() )
  
  # clean up
  system('ls sub*'); system('rm sub*')
  system('rm white_matter_mask.nii.gz')
  system('rm gray_matter_mask_out.nii.gz')
  system('rm white_matter_mask_or.nii.gz')
  system('rm gray_matter_mask_out_or.nii.gz')
  system('rm -R greyMatter')
  system('rm -R whiteMatter')
  system('rm -R volumetric')
  system('rm volumetricData*')
  system('rm whiteMatterLevelset.nii.gz')
  system('rm greyMatterLevelset.nii.gz')
  system('rm profileData_*')
  
  #move 2 folders up
  setwd('../..'); print( dir() )
  
  # copy data in SUMA folder
  filesToCopy <- dir( pattern='sub*' )
  for ( fIndx in 1:length( filesToCopy ) ) {
    instr <- sprintf( 'cp %s AHEAD_test/SUMA/%s', filesToCopy[fIndx], filesToCopy[fIndx] ); system( instr )
  }
  
  setwd( 'AHEAD_test' )
  setwd( 'SUMA' )
  print( dir() )
  
  # wm mask #linear resample and cutoff
  instr <- sprintf( '3dcalc -a aparc+aseg_REN_all.nii.gz -expr \u027within(a,0.95,1.05)+within(a,20.95,21.05)\u027 -prefix white_matter_mask_or.nii.gz'  ); system( instr )
  instr <- sprintf('3dresample -master %s -rmode Lin -inset white_matter_mask_or.nii.gz -prefix white_matter_mask.nii.gz', filesToCopy[1] ); system( instr )
  # gm mask 
  instr <- sprintf( '3dcalc -a aparc+aseg_REN_all.nii.gz -b lh.ribbon.nii -c rh.ribbon.nii -expr \u027within(a,0.95,1.05)+within(a,20.95,21.05)+b+c\u027 -prefix gray_matter_mask_out_or.nii.gz'  ); system( instr )
  instr <- sprintf('3dresample -master %s -rmode Lin -inset gray_matter_mask_out_or.nii.gz -prefix gray_matter_mask_out.nii.gz', filesToCopy[1] ); system( instr )  

  # volumetric call
  system( volumetricCodeCall )

  # t1 profile
  system( sprintf( 'cp %s inputData.nii.gz', filesToCopy[5] ) )
  # profile call
  system( profileCodeCall )
  # rename
  system('mv profilesOut_profiles.nii.gz profileData_t1w.nii.gz') #!!!
  # clean up
  system('rm inputData.nii.gz')
  
  # qsm profile
  system( sprintf( 'cp %s inputData.nii.gz', filesToCopy[1] ) )
  # profile call
  system( profileCodeCall )
  # rename
  system('mv profilesOut_profiles.nii.gz profileData_qsm.nii.gz')
  # clean up
  system('rm inputData.nii.gz')

  # r1m profile
  system( sprintf( 'cp %s inputData.nii.gz', filesToCopy[2] ) )
  # profile call
  system( profileCodeCall )
  # rename
  system('mv profilesOut_profiles.nii.gz profileData_r1m.nii.gz')
  # clean up
  system('rm inputData.nii.gz')

  # r2star profile
  system( sprintf( 'cp %s inputData.nii.gz', filesToCopy[3] ) )
  # profile call
  system( profileCodeCall )
  # rename
  system('mv profilesOut_profiles.nii.gz profileData_r2star.nii.gz')
  # clean up
  system('rm inputData.nii.gz')

  # t1map profile
  system( sprintf( 'cp %s inputData.nii.gz', filesToCopy[4] ) )
  # profile call
  system( profileCodeCall )
  # rename
  system('mv profilesOut_profiles.nii.gz profileData_t1Map.nii.gz')
  # clean up
  system('rm inputData.nii.gz')
  
  # t2star profile
  system( sprintf( 'cp %s inputData.nii.gz', filesToCopy[6] ) )
  # profile call
  system( profileCodeCall )
  # rename
  system('mv profilesOut_profiles.nii.gz profileData_t2Star.nii.gz')
  # clean up
  system('rm inputData.nii.gz')
  
}