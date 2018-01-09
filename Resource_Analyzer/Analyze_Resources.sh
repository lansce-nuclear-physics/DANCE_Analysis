#!/bin/bash
#set -x

THISHOST=$(hostname)

xterm -hold -e python Resource_Analyzer.py iowait_$THISHOST.txt &

#Start the iowait watcher
iostat -x 1 | tee iowait_$THISHOST.txt

#xterm -hold -e python Resource_Analyzer.py iowait_$THISHOST.txt

#Start the root plotter
#xterm -hold -e root Resource_Analyzer.C 
