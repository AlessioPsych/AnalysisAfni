#!/bin/bash

ANATFOLDER=$1

cd /analyse/Project0226/KastnerModel/anatomies_KastnerClassic
cd $ANATFOLDER

suma -niml &
sleep 2;

DriveSuma -echo_edu -com show_surf -label ICO -i_1D surfaces_folder_right/inflated_boundary01.1D.coord surfaces_folder_right/inflated_01_or.1D.topo  -surf_winding ccw

DriveSuma -echo_edu -com surf_cont -load_dset zzz_phase_encoded_wAverage_10_correct_anat_add_interp_surfaces_folder_right/boundary01_sm_zzz_phase_encoded_wAverage_10_correct_anat_add_surf.1D.dset -surf_label ICO -view_surf_cont y

#DriveSuma -echo_edu -com surf_cont -switch_dset Convexity10#DriveSuma -echo_edu -com surf_cont T_val 5

DriveSuma -echo_edu -com surf_cont -switch_dset boundary01_sm_zzz_phase_encoded_wAverage_10_correct_anat_add_surf.1D.dset

DriveSuma -echo_edu -com surf_cont -I_sb 8 -T_sb 7 -T_val 0.35

DriveSuma -echo_edu -com surf_cont -Clst -1 25 -UseClst y

DriveSuma -com surf_cont -switch_cmap Color_circle_ZSS

DriveSuma -echo_edu -com surf_cont -I_range 0.1 6.28 

cd /analyse/Project0226/KastnerModel/code

