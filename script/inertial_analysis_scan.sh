#!/bin/bash
files=$(find . -name "*_inertial.txt")
for i in ${files[@]}; do
	ii=$(basename $i)
	echo "--------- Inertial Analysis : $ii"
	inertial_stats.exe $ii
	echo
done
