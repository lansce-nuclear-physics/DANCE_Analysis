#!/bin/bash

#PATHTODATA=caen2018
#PATHTODATA=stage0_bin_automated
#PATHTODATA=stage0_bin
#PATHTODATA=dancedaq7
#PATHTODATA=/mnt/dance-data/2/caen2018
PATHTODATA=/mnt/dance-data/6/caen2018

#for i in `seq 67339 67357`;  #Co59 141.0mg
#for i in `seq 67359 67835`;  #Co59 329.4mg
#for i in `seq 68128 68900`;  #Co59 329.4mg
#for i in `seq 67836 68039`;  #Co59 9.9mg

#for i in `seq 60834 61251`;  #Cu65
#for i in `seq 63015 63589`;  #Cu65

#for i in `seq 63015 63230`;  #Cu65
#for i in `seq 63231 63400`;  #Cu65
#for i in `seq 63401 63589`;  #Cu65

#for i in `seq 64137 64692`;  #Cu65
#for i in `seq 64989 65265`;  #Cu65

#for i in `seq 64693 64833`;  #Cu63

#for i in `seq 60409 60451`;  #Au197 4mm
#for i in `seq 5010 5240`;  #Au197 4mm
#for i in `seq 5241 5410`;  #Au197 4mm
#for i in `seq 5410 5500`;  #Au197 4mm

#for i in `seq 5501 5600`;  #Pb208 
#for i in `seq 5601 5820`;  #Pb208 
#for i in `seq 5821 5940`;  #Pb208 

#for i in `seq 5841 6020`;  #Pb208 
#for i in `seq 6021 6160`;  #Pb208 

#for i in `seq 7001 7090`;  #Pb208 
#for i in `seq 7091 7440`;  #Pb208 

#for i in `seq 64834 64987`;  #Pb208 
#for i in `seq 68040 68103`;  #Pb208 

#for i in `seq 67519 67835`; 
#for i in `seq 67966 68039`;
#for i in `seq 68247 68900`;

for i in `seq $1 $2`;
#list="thresholdRuns.txt"
#for i in $(cat $list);

do
    echo Processing Run: $i;
  # ./DANCE_Analysis $PATHTODATA $i cfg_files/stage1_Pt194.cfg;
   ./DANCE_Analysis $PATHTODATA $i stage0_caen2018.cfg;

done
