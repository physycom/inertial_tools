#!/bin/bash

if (( $# < 1 )) 
then 
 echo "Usage: $0 filename "
 exit
fi

bin_fraction=0.01

{
  echo "## Shift scan @ $bin_fraction"
  echo "# shift # Entropy AX # Entropy AY # Entropy GZ # Mutual AY-AX # AY-GZ #"
} > scan_shift.dat

shift_values=( 10 1000 )
shift_values=$(seq 0 100 20000)

for shift in ${shift_values[@]}; do
  mutual_entropy.exe -bf $bin_fraction -s $shift $1 > scan_shift.log
  cat scan_shift.log
    
  values=($(cat scan_shift.log | grep = ))
  
  echo -e "$shift\t${values[3]}\t${values[7]}\t${values[11]}\t${values[15]}\t${values[19]}\t${values[23]}" >> ${1}_scan_shift.dat
done
