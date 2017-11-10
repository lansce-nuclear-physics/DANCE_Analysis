#!/bin/bash

#These all live on data/14/
for i in `seq 64693 64833`;  #Cu63

do
    echo Processing Run: $i;
    ./DANCE_Analysis $i
done
