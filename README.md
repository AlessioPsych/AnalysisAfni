# AnalysisAfni

# for nighresaddon on the glasgow server: run_nighres_command.sh from python3 to python / segmentScriptinhouse_dec... line 42: cp cruise_whole/test_cruise_whole_cruise-cortex.nii.gz cruise_cortex.nii.gz / move surfacesHemisphere.sh from layout to coregistration and change the atlas definition at the beginning



to set up:

to setup you need to install afni and R on an ubuntu machine (https://afni.nimh.nih.gov/pub/dist/doc/htmldoc/background_install/install_instructs/index.html)

once installed you need to download my toolbox (https://github.com/AlessioPsych/AnalysisAfni) into a folder (generally  in home/Programs), using the command:

git clone https://github.com/AlessioPsych/AnalysisAfni.git

you need to copy the file startAfniToobox_example01 from the toolbox directory

and change the following lines, with the path to your AnalysisAfni folder:

line 3: export AFNI_TOOLBOXDIR=/home/alessiof/Programs/analysisAfni # change to your folder
 line 4: PATH=$PATH:/home/alessiof/Programs/analysisAfni # change to your folder

and change the following lines, with the path to your afni installation folder:

line 23 export AFNI_INSTALLDIR=/usr/local/bin/ # change to AFNI install folder (bin)
 line 23 export AFNI_ATLASDIR=/usr/local/bin/ # change to AFNI folder where the atlases live


you can find the afni installation folder typing on the terminal:

locate 3dSkullStrip

