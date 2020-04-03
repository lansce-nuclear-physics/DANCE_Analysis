
prefix="Stage1_Histograms_Run_"
dir="stage1_root/"
prefout="_Stage1_Histograms_Runs_"
suffix="_5ns_CW_0ns_CBT_0ns_DEBT.root "
inter="_"
for i in $(seq $1 $2)
do
  list+=$dir$prefix${i}*$suffix
done 

hadd -k -f $dir$3$prefout$1$inter$2$suffix $list
