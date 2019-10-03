#!/bin/bash
if [ "$#" -ne 3 ]; then
  echo "missing command line argument"
  exit
fi
for (( i=$2; i<$3; i++))
do
  # echo Processing Run: $i from $1"
  ./DANCE_Analysis $1 $i stage0_caen2018.cfg;
    #./DANCE_Analysis $1 $i cfg_files/stage1_Pt194.cfg
     #./DANCE_Analysis $1 $i stage0_caen2015.cfg;
   # ./DANCE_Analysis $1 $i cfg_files/stage1_Au197.cfg
   # ./DANCE_Analysis $1 $i stage1_Eu152.cfg
#   ./DANCE_Analysis $1 $i cfg_files/stage1_Co59.cfg
  # ./DANCE_Analysis $1 $i cfg_files/stage1_Cu65.cfg
done
