#!/bin/bash

for (( i=$2; i<$3; i++))
do
  # echo Processing Run: $i;
   # echo  $1"  "$i"  stage0.cfg";
    ./DANCE_Analysis $1 $i stage0.cfg;
    ./DANCE_Analysis $1 $i stage1.cfg
done