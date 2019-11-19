#assumes you're running this in the DANCE_Analysis directory, where the script nominally lives...

R0="stage0_root"
B0="stage0_bin"
R1="stage1_root"
B1="stage1_bin"
S0="stage0_simulated"

for DIRECTORY in {$R0,$B0,$R1,$B1,$S0}; do
  echo "setting up" $DIRECTORY
  if [ ! -d $DIRECTORY ]; then
    mkdir $DIRECTORY
  fi
done

