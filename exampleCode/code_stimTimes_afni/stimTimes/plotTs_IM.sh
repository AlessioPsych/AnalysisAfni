3dDeconvolve -nodata 350 1 -polort -1 -num_stimts 1 \
		 -basis_normall 1 \
                 -stim_times_IM 1 q_AM_short.1D 'dmUBLOCK(1)'      \
  -global_times                   \
                 -x1D stdout: |                         \
     1dplot -stdin -thick                               
