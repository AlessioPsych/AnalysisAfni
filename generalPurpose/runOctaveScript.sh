#!/usr/bin/octave -qf
# example octave script

printf ("%s", program_name ());
arg_list = argv ();
addpath ( genpath ( arg_list{1} ) )
for i = 1:nargin
  printf (" %s", arg_list{i});
endfor

a = BrikLoad( 'DEPTH_al_epi.nii.gz' )
printf ("\n");

# printf("Elapsed time: %.4f seconds", elapsed);
