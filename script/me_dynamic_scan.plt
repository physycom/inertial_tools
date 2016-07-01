#!/gnuplot
FILE_IN='scan_dyn.dat'
FILE_OUT='scan_dyn.png'
set terminal pngcairo size 1280,720 font ",20"
set output FILE_OUT
set title 'Mutual Entropy'
set xlabel 'bin record number'
set x2label 'bin fraction'
set ylabel 'Entropy'
set y2label 'ratio'
set ytics nomirror
set y2tics
set x2tics
plot FILE_IN u ($1*$7):5 w lines lt 1 lc rgb 'blue' lw 3 t 'MI(a_x , g_z)' axes x1y1,\
     FILE_IN u ($1*$7):6 w lines lt 1 lc rgb 'red' lw 3 t 'MI(a_y , g_z)' axes x1y1,\
     FILE_IN u 1:($5/$6) w lines lt 1 lc rgb 'green' lw 3 t 'ratio (x/y)' axes x2y2

