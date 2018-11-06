import os
import sys
import nighres

#example call
# python delMe.py aaa bbb


#out_dir = os.path.join(os.getcwd(), '/home/fracasso/programs/nighres/delMe')

#dataset = nighres.data.download_7T_TRT(out_dir)

print "This is the name of the script: ", sys.argv[0]
print "Number of arguments: ", len(sys.argv)
print "The arguments are: " , str(sys.argv)

nighres.surface.probability_to_levelset(sys.argv[1], save_data=sys.argv[2], output_dir=sys.argv[3], file_name=sys.argv[4])
