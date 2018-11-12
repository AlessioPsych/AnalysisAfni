args <- commandArgs(T)
print( args )

setwd('/home/fracasso/data/tests/MP2RAGE_test_nighres')

args <- c('contrast_image1=MPRAGE_ss.nii', 'contrast_type1=Mp2rage7T',
          'contrast_image2=MPRAGE_T1_ss.nii', 'contrast_type2=T1map7T', 
          'save_data=True', 'file_name=subj_sess', 'output_dir=/home/fracasso/data/tests/MP2RAGE_test_nighres')

mainDir <- getwd()

nArgs <- length(args)
for (k in 1:nArgs) {
  str
  if (k==1) { instr <- sprintf('nighres.brain.mgdm_segmentation(%s, ', args[k]) }
  if (k>1 & k<nArgs) {
    instr <- sprintf( '%s %s, ', instr, args[k] )
  } 
  if (k==nArgs) {
    instr <- sprintf( '%s %s)', instr, args[k] )
  }
}
instr

instr <- sprintf( 'import' )


#mgdm_results = nighres.brain.mgdm_segmentation(
#  contrast_image1=skullstripping_results['t1w_masked'],
#  contrast_type1="Mp2rage7T",
#  contrast_image2=skullstripping_results['t1map_masked'],
#  contrast_type2="T1map7T",
#  save_data=True, file_name="sub001_sess1",
#  output_dir=out_dir)

