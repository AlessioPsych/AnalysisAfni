args <- commandArgs(T)
print( args )

#args <- c( 'A1_T123DEPI.nii', 'A1_T123DEPITU.nii', 'T123DEPI_unif_mask_dil.nii', 'T123DEPITU_unif_mask_dil.nii', '5' )
#setwd('/data2/scdropbox/Students/Tessa/DATA/Participants/A_Imagery/20170214_A1_Cerebellum/Anatomical/T123DEPI_applyTU')

# get actual dir
mainDir <- getwd()

t1 <- args[1]
t1TU <- args[2]
#secondInv <- args[3]
#secondInvTU <- args[4]

#instr <- sprintf( '3dAutomask -clfrac 0.6 -prefix secondInv_mask.nii.gz %s', secondInv );
#system( instr )
#instr <- sprintf( '3dAutomask -clfrac 0.6 -prefix secondInv_mask_TU.nii.gz %s', secondInvTU );
#system( instr )

instr <- sprintf( '3dcalc -a %s -b %s -expr \u0027a*step(b)\u0027 -prefix t1_mask_no_blur.nii.gz', t1, args[3] );
system( instr )
instr <- sprintf( '3dcalc -a %s -b %s -expr \u0027a*step(b)\u0027 -prefix t1_TU_mask_no_blur.nii.gz', t1TU, args[4] );
system( instr )

instr <- sprintf('3dQwarp -source t1_mask_no_blur.nii.gz -base t1_TU_mask_no_blur.nii.gz -prefix warpTop -verb -iwarp -pblur 0.05 0.05 -blur -1 -1 -noweight -minpatch %s -plusminus', args[5] )
system( instr )

# storage directory
print('###################')
print('###################')
print('clean up transform')
print('###################')
print('###################')
targetDir <- 'topUpDir_T1w'
flagDir <- dir.create( file.path(mainDir, targetDir) )
if (flagDir==FALSE) { # directory already exists!
  msg <- sprintf( 'Remove the directory %s_folder to proceed', targetDir )
  warning( msg )
  stopifnot(flagDir)  
}

#system('rm secondInv_mask.nii.gz')
#system('rm secondInv_mask_TU.nii.gz')
system('mv t1_mask_no_blur.nii.gz topUpDir_T1w/')
system('mv t1_TU_mask_no_blur.nii.gz topUpDir_T1w/')
system('mv warpTop* topUpDir_T1w/')
#instr <- sprintf('3dQwarp -source t1_mask_no_blur.nii.gz -base t1_TU_mask_no_blur.nii.gz -prefix SUBJ_TU_no_blur -verb -iwarp -blur 0 0 -plusminus');
#system( instr )

#instr = sprintf('3dNwarpApply -master t1_mask_no_blur.nii.gz -source t1_mask_no_blur.nii.gz -nwarp SUBJ_TU_no_blur_PLUS_WARP+orig -interp NN -prefix dataWarped_no_blur+orig');
#system( instr )



#instr = sprintf('3dcalc -a secondInv_mask.nii.gz -b secondInv_mask_TU.nii.gz -expr ''step((a+b)-1.1)*1000'' -prefix maskVolume.nii.gz');
#system( instr )
#instr = sprintf('3dmerge -1blur_sigma 4.0 -prefix maskVolume_blur.nii.gz maskVolume.nii.gz');
#system( instr )
#instr = sprintf('3dcalc -a maskVolume_blur.nii.gz -expr ''step(a-100)'' -prefix maskVolume_blur_thr.nii.gz');
#system( instr )

#instr = sprintf( '3dcalc -a %s -b maskVolume_blur_thr.nii.gz -expr ''a*step(b)'' -prefix t1_mask.nii.gz', t1 );
#system( instr )
#instr = sprintf( '3dcalc -a %s -b maskVolume_blur_thr.nii.gz -expr ''a*step(b)'' -prefix t1_TU_mask.nii.gz', t1TU );
#system( instr )

#instr = sprintf('3dQwarp -source t1_mask.nii.gz -base t1_TU_mask.nii.gz -prefix SUBJ_TU -verb -iwarp -blur 0 0 -plusminus');
#system( instr )

#instr = sprintf('3dNwarpApply -master t1_mask.nii.gz -source t1_mask.nii.gz -nwarp SUBJ_TU_PLUS_WARP+orig -interp NN -prefix dataWarped+orig');
#system( instr )
