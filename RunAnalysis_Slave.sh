#!/bin/bash

for (( i=$2; i<$3; i++))
do
  # echo Processing Run: $i from $1"
   # ./DANCE_Analysis $1 $i stage0_caen2018.cfg;
    #./DANCE_Analysis $1 $i stage0_caen2015.cfg;
  #  ./DANCE_Analysis $1 $i stage1_Er167.cfg
   # ./DANCE_Analysis $1 $i stage1_Eu152.cfg
    ./DANCE_Analysis $1 $i stage1_Cd114.cfg
done