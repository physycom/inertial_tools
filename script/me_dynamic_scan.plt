#!/gnuplot
FILE_IN='scan.dat'
FILE_OUT='scan.png'
set terminal pngcairo size 1280,720 font ",20"
set output FILE_OUT
set title 'Mutual Entropy'
set xlabel 'bin fraction'
set ylabel 'Entropy'
set y2label 'ratio'
set ytics nomirror
set y2tics
plot FILE_IN u 1:5 w lines lt 1 lc rgb 'blue' lw 3 t 'MI(a_y , a_x)' axes x1y1,\
     FILE_IN u 1:6 w lines lt 1 lc rgb 'red' lw 3 t 'MI(a_y , g_z)' axes x1y1,\
     FILE_IN u 1:($5/$6) w lines lt 1 lc rgb 'green' lw 3 t 'ratio' axes x1y2

