#!/bin/bash
#set -x

rsync -av /mnt/hygelac-data/29/dance/Automated_DANCE_Analysis_Results/TimeDeviations/TimeDeviations_Run_10{4,5,6,8}*.txt ../TimeDeviations/

rsync -av /home/cfry/DANCE_Alpha_Calibrator/ParamOutput ../Calibrations/
