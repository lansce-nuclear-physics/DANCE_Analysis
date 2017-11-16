#!/bin/bash

PATHTODATA=MIDAS_Data
#PATHTODATA=MIDAS_Data_dancedaq_6

#These all live on data/14/
#for i in `seq 64693 64833`;  #Cu63

#for i in `seq 63015 63589`;  #Cu65

#for i in `seq 64989 65265`;  #Cu65
#for i in `seq 65047 65265`;  #Cu65

#for i in `seq 64137 64692`;  #Cu65

#for i in `seq 64834 64987`;  #Pb208



#for i in `seq 67339 67357`;  #Co59 141.0mg
for i in `seq 67366 67401`;  #Co59 329.4mg


do
    echo Processing Run: $i;
    ./DANCE_Analysis $PATHTODATA $i;
done
