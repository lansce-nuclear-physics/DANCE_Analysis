#!/bin/bash
#set -x

rsync -av /mnt/hygelac-data/29/dance/Automated_DANCE_Analysis_Results/TimeDeviations/TimeDeviations_Run_108*.txt ../TimeDeviations/

rsync -av /mnt/hygelac-data/29/dance/Automated_DANCE_Analysis_Results/Calibrations/param*108*.txt ../Calibrations/