unlink stage0_root
unlink stage0_bin
unlink stage0_simulated
unlink stage0_root_automated1
unlink stage0_bin_automated1
unlink stage0_root_automated2
unlink stage0_bin_automated2
unlink stage1_root
unlink stage1_bin

for DIRECTORY in data*/; do
  unlink ${DIRECTORY%\/}
done

for DIRECTORY in hygelac*/; do
  unlink ${DIRECTORY%\/}
done
