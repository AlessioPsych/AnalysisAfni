#!/bin/bash

ANATFOLDER=$1

cd /analyse/Project0226/KastnerModel/anatomies_KastnerModel
cd $ANATFOLDER

suma -niml &
sleep 2;

DriveSuma -echo_edu -com show_surf -label ICO -i_1D surfaces_folder_left/inflated_boundary01.1D.coord surfaces_folder_left/inflated_01_or.1D.topo  -surf_winding ccw

DriveSuma -echo_edu -com surf_cont -load_dset zzz_kastner_model_params_anat_add_interp_surfaces_folder_left/boundary01_sm_zzz_kastner_model_params_anat_add_surf.1D.dset -surf_label ICO -view_surf_cont y

DriveSuma -echo_edu -com surf_cont -switch_dset boundary01_sm_zzz_kastner_model_params_anat_add_surf.1D.dset

DriveSuma -echo_edu -com surf_cont -I_sb 11 -T_sb 14 -T_val 0.15

DriveSuma -echo_edu -com surf_cont -Clst -1 25 -UseClst y

DriveSuma -com surf_cont -switch_cmap Color_circle_ZSS

DriveSuma -echo_edu -com surf_cont -I_range 0.1 6.28 

#DriveSuma -echo_edu -com surf_cont -switch_dset Convexity

cd /analyse/Project0226/KastnerModel/code

