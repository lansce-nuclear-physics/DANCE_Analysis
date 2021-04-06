#!/bin/bash
#set -x

#./Copy_Automated.sh

#where the data are
input=hygelac13
#input=./stage0_bin_automated
#input=./stage0_bin
#sh gates2019.sh
#input=data6
#input=data4

#Get the start and end run from the command line args
startrun=$1;
endrun=$2;

#number of xterms with parallel jobs running
nxterms=$3   #number of xterms to launch at any given time

#make sure there are enough command line arguments
if [ "$#" -ne 3 ]; then
  echo "missing command line argument (startrun, endrun, nxterms)"
#  exit
fi

#make sure the start run is less than the end run
if [ "$endrun" -lt "$startrun" ]
then
    echo "Problem with start and end run values"
    exit
fi

#the start and end are inclusive so there is a +1
nruns=$(expr $endrun - $startrun + 1)

echo "Analyzing: "$nruns" Runs" 

nrunsperxterm=$(expr $nruns / $nxterms ) #+ 1) need the +1 if very assymetric

echo "Starting: "$nxterms" terminals with "$nrunsperxterm" runs each"

#Start the various analysis jobs
counter=0;
currentstart=$startrun

for (( j=0; j<$(expr $nxterms - 1); j++))
do
    currentstart=$(( $startrun + $(expr $counter \* $nrunsperxterm ) ))
    counter=$(( $counter+1 ))  
    currentstop=$(expr $startrun + $(expr $counter \* $nrunsperxterm))
   # echo $currentstart"  "$currentstop

    #xterm -e ./RunAnalysis_Minion.sh $input $currentstart $currentstop &
     ./RunAnalysis_Minion.sh $input $currentstart $currentstop $nucleus &


    sleep 10
done

#the last one gets the remainder
currentstart=$(( $startrun + $(expr $counter \* $nrunsperxterm ) ))
currentstop=$(( $endrun + 1 ))
# echo $currentstart"  "$currentstop

#xterm -e ./RunAnalysis_Minion.sh $input $currentstart $currentstop &
./RunAnalysis_Minion.sh $input $currentstart $currentstop $nucleus &


exit 0
