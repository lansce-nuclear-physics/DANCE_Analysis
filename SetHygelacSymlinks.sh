#assumes you're running this in the DANCE_Analysis directory, where the script nominally lives...
#if you have directories on multiple volumes and want to specify, pass it as command line arg

UDIR=${HOME#/home/}
UDIR=$UDIR"/"
HYGELAC="/mnt/hygelac-data/"
DANCE="/dance/"

ISDIR=0

for VOLUME in "20" "10"; do
  if [ -d $HYGELAC$VOLUME$DANCE$UDIR ]; then
    echo "setting up things in volume" $VOLUME
    ISDIR=+1
    break
  fi
done

if [ $ISDIR -eq 0 ]; then
  mkdir $HYGELAC$VOLUME$DANCE$UDIR
fi

if [ $# -ge 1 ]; then
  VOLUME=$1
fi

if [ ! -d $HYGELAC$VOLUME$DANCE$UDIR ]; then
  echo "talk to someone about getting added to the dance group"
  exit 1
fi
 
ln -s /mnt/hygelac-data/31/dance/Automated_DANCE_Analysis_Results/stage0_bin stage0_bin_automated
ln -s /mnt/hygelac-data/31/dance/Automated_DANCE_Analysis_Results/stage0_root stage0_root_automated

R0="stage0_root"
B0="stage0_bin"
R1="stage1_root"
B1="stage1_bin"
S0="stage0_simulated"

for DIRECTORY in  $R0 $B0 $R1 $B1 $S0 ; do
  echo "setting up" $DIRECTORY
  if [ ! -d $HYGELAC$VOLUME$DANCE$UDIR$DIRECTORY ]; then
    mkdir $HYGELAC$VOLUME$DANCE$UDIR$DIRECTORY
  fi
  ln -s $HYGELAC$VOLUME$DANCE$UDIR$DIRECTORY $DIRECTORY
  if [ ! -d $DIRECTORY ]; then
    echo "something has gone wrong, ask people about permissions"
  fi
done

echo "setting up links to dancedaq data"
DDAQ="/mnt/dance-data/"
DDEND="/caen2018"
for DVOL in $(seq 1 7); do 
  ln -s $DDAQ$DVOL$DDEND data$DVOL 
done

echo "setting up links to hygelac data"
for HVOL in $(seq 11 14); do
  ln -s $HYGELAC$HVOL$DANCE$DDEND hygelac$HVOL
done

