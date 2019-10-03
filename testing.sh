for i in $(cat ./runlist.txt); do
  ./DANCE_Analysis hygelac24 $i stage0_caen2018.cfg 
  echo $i
  read -n 1 -s -r -p "last run was"
done
