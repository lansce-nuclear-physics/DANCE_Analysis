#!/bin/bash

source ~/.bashrc

(
# Uncomment the following line for debugging output
#set -x

/bin/date +%T;

umask -S
umask 022
umask -S

#FIRST_RUN=67339
#LAST_RUN=68115

#FIRST_RUN=60622
#LAST_RUN=61251

#FIRST_RUN=62695
#LAST_RUN=62735

#FIRST_RUN=63015
#LAST_RUN=65265

#FIRST_RUN=59843
#LAST_RUN=60584

#These are on hygelac13
FIRST_RUN=59236
LAST_RUN=59679

#FIRST_RUN=59680
#LAST_RUN=59842

#Cu63
#FIRST_RUN=64693
#LAST_RUN=64833

#Cu65
#FIRST_RUN=63015
#LAST_RUN=63589

#Pb208
#FIRST_RUN=60622
#LAST_RUN=60811

#Au197
#FIRST_RUN=60409
#LAST_RUN=60451

#Au197
#FIRST_RUN=60572
#LAST_RUN=60584

#Au197
#FIRST_RUN=63590
#LAST_RUN=63649

#Au197
#FIRST_RUN=60409
#LAST_RUN=60451

DATA_DIR=/mnt/hygelac-data/13/dance/caen2015/
ANALYSIS_DIR=/home/cprokop/CJP/DANCE_Analysis/
STAGE0_DIR=$ANALYSIS_DIR/stage0_root/
STAGE1_DIR=$ANALYSIS_DIR/stage1_root/
ALPHA_CALIB_DIR=/home/cprokop/CJP/DANCE_Alpha_Calibrator

for FILE in $DATA_DIR/run*.mid.gz ; do
    
    if [ -s $FILE ] ; then 
	GZ_FILE=`basename $FILE`
	BASENAME=${GZ_FILE%.mid.gz}
	RUN=${BASENAME#run}
	
        # This is weird--to deal with the fact that there are leading zeros,
        #   we pretend the numbers are in hex so that bash can do a numeric
        #   comparison
	HEXRUN=0x$RUN
	HEXFIRST_RUN=0x$FIRST_RUN
	HEXLAST_RUN=0x$LAST_RUN
	if (( $HEXRUN > $HEXLAST_RUN )) ; then
	    echo "Last run " $LAST_RUN "completed"
	    exit 0
	fi

	if (( $HEXRUN >= $HEXFIRST_RUN )) ; then 
	    #check for existing rootfile
	    if [ ! -f $STAGE0_DIR/Stage0_Histograms_Run_${RUN#0}.root ] ; then
		#switch to the analysis directory 
		cd $ANALYSIS_DIR
		
		#Check to see if fare is built and then run the analysis
		if [ -f DANCE_Analysis ] ; then
		    echo "about to run stage 0 analyis"
		./DANCE_Analysis $DATA_DIR $RUN stage0.cfg

		fi		
	    fi
	    
	    #switch to alpha calibrations
	    if [ ! -f $ALPHA_CALIB_DIR/ParamOutput/param_out_${RUN#0}.txt ] ; then
		cd $ALPHA_CALIB_DIR
		if [ -f DANCE_Alpha_Calibrator ] ; then
		    if [ -f $STAGE0_DIR/Stage0_Histograms_Run_${RUN#0}.root ] ; then
			echo "about to run calibrations"
			./DANCE_Alpha_Calibrator $RUN
			cp $ALPHA_CALIB_DIR/ParamOutput/param_out_${RUN#0}.txt $ANALYSIS_DIR/Calibrations/.
		    fi
		fi
	    fi   

	    #check for existing rootfile
	    if [ ! -f $STAGE1_DIR/Stage1_Histograms_Run_${RUN#0}.root ] ; then
		#switch to the analysis directory 
		cd $ANALYSIS_DIR
		
		#Check to see if fare is built and then run the analysis
		if [ -f DANCE_Analysis ] ; then
		    echo "about to run stage 1 analyis"
		./DANCE_Analysis $DATA_DIR $RUN stage1.cfg

		fi		
	    fi
	fi
    fi
done

) | tee -a /tmp/`basename $0`.`date +%F`.log


