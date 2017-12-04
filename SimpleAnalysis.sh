#!/bin/bash

PATHTODATA=hygelac14
#PATHTODATA=MIDAS_Data_dancedaq_6

#These all live on data/14/
#for i in `seq 64693 64833`;  #Cu63

#for i in `seq 63015 63589`;  #Cu65

#for i in `seq 64989 65265`;  #Cu65
#for i in `seq 65047 65265`;  #Cu65

#for i in `seq 64137 64692`;  #Cu65

#for i in `seq 64834 64987`;  #Pb208



#for i in `seq 67339 67357`;  #Co59 141.0mg
#for i in `seq 67359 67399`;  #Co59 329.4mg
#for i in `seq 67400 67449`;  #Co59 329.4mg
#for i in `seq 67450 67499`;  #Co59 329.4mg
#for i in `seq 67500 67549`;  #Co59 329.4mg
#for i in `seq 67550 67599`;  #Co59 329.4mg
#for i in `seq 67600 67649`;  #Co59 329.4mg
#for i in `seq 67650 67699`;  #Co59 329.4mg
#for i in `seq 67700 67760`;  #Co59 329.4mg

#for i in `seq 60409 60451`;  #Au197 4mm
for i in `seq 60572 60584`;  #Au197 4mm
#for i in `seq 63590 63649`;  #Au197 4mm
#for i in `seq 68104 68126`;  #Au197 4mm


do
    echo Processing Run: $i;
    ./DANCE_Analysis $PATHTODATA $i stage1.cfg;
done