#!/usr/bin/gnuplot

set terminal pngcairo
set output "PERBX_fsgrid.png"

set xlabel "y (cells)"
set ylabel "z (cells)"

set xrange [0:5]
set yrange [0:5]
set cbrange [-1:1]

plot \
   "PERBX_fsgrid.dat" using ($1+.5):($2+.5):3 matrix with image title "PERBX",\
   'samples.dat' using 2:3:4 lc palette title "Samples"
