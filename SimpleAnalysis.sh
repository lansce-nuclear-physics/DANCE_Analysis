#!/bin/bash

PATHTODATA=hygelac14
#PATHTODATA=MIDAS_Data_dancedaq_6

#These all live on data/14/

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
#for i in `seq 60572 60584`;  #Au197 4mm
#for i in `seq 63590 63649`;  #Au197 4mm
#for i in `seq 68104 68126`;  #Au197 4mm

#for i in `seq 59715 59842`;  #Pb208 
#for i in `seq 60622 60811`;  #Pb208 
#for i in `seq 62695 62735`;  #Pb208 
#for i in `seq 63650 64136`;  #Pb208 
#for i in `seq 64834 64987`;  #Pb208 
for i in `seq 68040 68103`;  #Pb208 

do
    echo Processing Run: $i;
    ./DANCE_Analysis $PATHTODATA $i stage0.cfg;
    ./DANCE_Analysis $PATHTODATA $i stage1.cfg;

done