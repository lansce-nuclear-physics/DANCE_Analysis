#!/bin/bash
#set -x

./Copy_Automated.sh

#where the data are
#input=./hygelac22/
input=./stage0_bin_automated
#input=./dancedaq3/
#input=./stage0_bin

#input=./dancedaq7/

#Get the start and end run from the command line args
startrun=$1;
endrun=$2;

#number of xterms with parallel jobs running
nxterms=$3   #number of xterms to launch at any given time

#make sure the start run is less than the end run
if [ "$endrun" -lt "$startrun" ]
then
    echo "Problem with start and end run values"
    exit
fi

#the start and end are inclusive so there is a +1
nruns=$(expr $endrun - $startrun + 1)

echo "Analyzing: "$nruns" Runs" 

nrunsperxterm=$(expr $nruns / $nxterms)

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

    #xterm -e ./RunAnalysis_Slave.sh $input $currentstart $currentstop &
     ./RunAnalysis_Slave.sh $input $currentstart $currentstop &


    sleep 10
done

#the last one gets the remainder
currentstart=$(( $startrun + $(expr $counter \* $nrunsperxterm ) ))
currentstop=$(( $endrun + 1 ))
# echo $currentstart"  "$currentstop

#xterm -e ./RunAnalysis_Slave.sh $input $currentstart $currentstop &
./RunAnalysis_Slave.sh $input $currentstart $currentstop &

exit 0
