#!/bin/bash

ANATFOLDER=$1

cd /analyse/Project0226/KastnerModel/anatomies_KastnerClassic
cd $ANATFOLDER

suma -niml &
sleep 2;

DriveSuma -echo_edu -com show_surf -label ICO -i_1D surfaces_folder_left/inflated_boundary01.1D.coord surfaces_folder_left/inflated_01_or.1D.topo  -surf_winding ccw

DriveSuma -echo_edu -com surf_cont -load_dset ppp_prf_nocss_params_anat_add_interp_surfaces_folder_left/boundary01_sm_ppp_prf_nocss_params_anat_add_surf.1D.dset -surf_label ICO -view_surf_cont y

DriveSuma -echo_edu -com viewer_cont -load_view viewRightParietal.niml.vvs

DriveSuma -echo_edu -com surf_cont -switch_dset Convexity

DriveSuma -echo_edu -com surf_cont -I_sb 0 -T_sb 0 -T_val 5

DriveSuma -echo_edu -com surf_cont -switch_dset boundary01_sm_ppp_prf_nocss_params_anat_add_surf.1D.dset

DriveSuma -echo_edu -com surf_cont -I_sb 22 -T_sb 17 -T_val 0.20

DriveSuma -echo_edu -com surf_cont -Clst -1 25 -UseClst y

DriveSuma -com surf_cont -load_cmap /analyse/Project0226/KastnerModel/code/bensonLeft.niml.cmap #Color_circle_ZSS

#DriveSuma -com surf_cont -load_cmap Color_circle_ZSS

DriveSuma -echo_edu -com surf_cont -I_range -3.14 3.14 -Dim 0.6 

cd /analyse/Project0226/KastnerModel/code

