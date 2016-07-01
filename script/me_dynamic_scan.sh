#!/bin/bash

{
  echo "## Dynamic bin scan"
  echo "# bin_fraction # Entropy AX # Entropy AY # Entropy GZ # Mutual AY-AX # AY-GZ # samples #"
} > scan_dyn.dat

#bf_values=( 0.20 0.11 )
bf_values=$(seq 0.01 0.005 0.35)

for bf in ${bf_values[@]}; do
  mutual_entropy.exe -bf $bf data.txt > scan.log
  cat scan.log
  echo
  
  values=($(cat scan.log | grep = ))
  
  echo -e "$bf\t${values[3]}\t${values[7]}\t${values[11]}\t${values[15]}\t${values[19]}\t${values[23]}" >> scan_dyn.dat
done
