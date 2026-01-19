3dDeconvolve -nodata 350 1 -polort -1 -num_stimts 3 \
                 -stim_times_AM1 1 q_short.1D 'dmUBLOCK'      \
                 -stim_times_AM1 2 q_short.1D 'dmUBLOCK(1)'   \
                 -stim_times_AM1 3 q_short.1D 'dmUBLOCK(-4)'  \
                 -x1D stdout: |                         \
     1dplot -stdin -thick                               \
            -ynames 'dmUBLOCK' 'dmUB(1)' 'dmUB(-4)'        
