import os
import sys

print "This is the name of the script: ", sys.argv[0]
print "Number of arguments: ", len(sys.argv)
print "The arguments are: " , str(sys.argv)

print sys.argv[len(sys.argv)-2]
mainDir = sys.argv[len(sys.argv)-2]
os.chdir( mainDir )

print sys.argv[len(sys.argv)-1]
nighresDir = sys.argv[len(sys.argv)-1]
os.chdir( nighresDir )
#import nighres
os.chdir( mainDir )


import sys
os.chdir( '/home/fracasso/programs/nighres/cbstools' )
sys.path.insert(0, "/home/fracasso/programs/nighres/cbstools")

os.chdir( '/home/fracasso/programs/nighres/nighres' )

sys.path.insert(0, "/home/fracasso/programs/nighres/nighres")

import nighres

mgdm_results = nighres.brain.mgdm_segmentation(
                        contrast_image1=skullstripping_results['t1w_masked'],
                        contrast_type1="Mp2rage7T",
                        contrast_image2=skullstripping_results['t1map_masked'],
                        contrast_type2="T1map7T",
                        save_data=True, file_name="sub001_sess1",
                        output_dir=out_dir)

#copy temporary file in nighresDir and execut from there

#nighres.surface.probability_to_levelset(contrast_image1=sys.argv[1], contrast_type1=sys.argv[2], contrast_image1=sys.argv[2], contrast_type2=sys.argv[3],)
