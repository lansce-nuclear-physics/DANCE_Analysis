#I make no guarantees this is complete, but an attempt was made.
#Feel free to tell me I'm wrong or add granularity
if [ -z $1 ]; then
  echo "give me a run number"
  exit 1
fi

RUNNUM=$1
echo "setting up calibrations"
cd ../Calibrations
unlink calib_ideal.dat
#2018
if [ $RUNNUM -le 111878 ]; then
  ln -s ../DANCE_Analysis/calib_ideal_2018.dat calib_ideal.dat
#2019 normal
elif [ $RUNNUM -le 117260 ]; then
  ln -s ../DANCE_Analysis/calib_ideal_2019.dat calib_ideal.dat
#2019 Tl
elif [ $RUNNUM -le 117640 ]; then
  ln -s ../DANCE_Analysis/calib_ideal_Tl204.dat calib_ideal.dat
#2020 normal
elif [ $RUNNUM -le 118734 ]; then 
  ln -s ../DANCE_Analysis/calib_ideal_2019.dat calib_ideal.dat
#2020 Tl
elif [ $RUNNUM -le 119056 ]; then  
  ln -s ../DANCE_Analysis/calib_ideal_Tl204.dat calib_ideal.dat
#2020 normal
elif [ $RUNNUM -le 119796 ]; then
  ln -s ../DANCE_Analysis/calib_ideal_2019.dat calib_ideal.dat
#2020 neuance
else 
  ln -s ../DANCE_Analysis/calib_ideal_2019.dat calib_ideal.dat
#y'all presumably want to put something in here that includes neuance
fi  

cd ../DANCE_Analysis/Gates
echo "setting up gates"
unlink Alpha.dat
unlink Gamma.dat
unlink Retrigger.dat

#2018
if [ $RUNNUM -le 111878 ]; then
  ln -s Alpha_109430.dat Alpha.dat
  ln -s Gamma_109430.dat Gamma.dat
  ln -s Retrigger_105138.dat Retrigger.dat
#2019 normal
elif [ $RUNNUM -le 117260 ]; then
  ln -s Alpha_112775.dat Alpha.dat
  ln -s Gamma_113158.dat Gamma.dat
  ln -s Retrigger_113341.dat Retrigger.dat
#2019 Tl
elif [ $RUNNUM -le 117640 ]; then
  ln -s Alpha_tl.dat Alpha.dat
  ln -s Gamma_113158.dat Gamma.dat
  ln -s Retrigger_Tl.dat Retrigger.dat
#2020 normal
elif [ $RUNNUM -le 118734 ]; then 
  ln -s Alpha_112775.dat Alpha.dat
  ln -s Gamma_113158.dat Gamma.dat
  ln -s Retrigger_113341.dat Retrigger.dat
#2020 Tl
elif [ $RUNNUM -le 119056 ]; then  
  ln -s Alpha_tl.dat Alpha.dat
  ln -s Gamma_113158.dat Gamma.dat
  ln -s Retrigger_Tl.dat Retrigger.dat
#2020 normal
elif [ $RUNNUM -le 119796 ]; then
  ln -s Alpha_112775.dat Alpha.dat
  ln -s Gamma_113158.dat Gamma.dat
  ln -s Retrigger_113341.dat Retrigger.dat
#2020 neuance
else
  ln -s Alpha_112775.dat Alpha.dat
  ln -s Gamma_113158.dat Gamma.dat
  ln -s Retrigger_113341.dat Retrigger.dat
fi  

echo "setting up config"
cd ../Config

unlink DanceMap.txt
unlink TMatrix.txt

if [ $RUNNUM -le 108028 ]; then
  ln -s DanceMap_23Nov15.txt DanceMap.txt
elif [ $RUNNUM -gt 119796 ]; then
  ln -s DanceMap_19Nov2020_neuance.txt DanceMap.txt
else
  ln -s DanceMap_7Nov18.txt DanceMap.txt
fi

if [ $RUNNUM -gt 119796 ]; then
  ln -s TMatrix_2020_neuance.txt TMatrix.txt
else
  ln -s TMatrix_2019.txt TMatrix.txt
fi

cd ..
