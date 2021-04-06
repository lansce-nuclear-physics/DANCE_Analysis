for i in `cat jumbled.txt` 
do 
  echo $i
  ./DANCE_Analysis hygelac11 $i cfg_files/stage1_mid.cfg
done
