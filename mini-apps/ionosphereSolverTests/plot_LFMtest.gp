#!/usr/bin/gnuplot -persist

set terminal pngcairo size 500,500
set output "LFMtest.png"

set xlabel "Latitude θ"
set ylabel "Φ(θ, φ) / sin(φ) (kV)"

set key bottom center

plot [45:90] [0:120]  \
      'lfmtest.dat' using ($1/3.14159*180):($4/1e3*1.05)  pt 0 title "our result", \
      'merkinReference.dat' using (90 - $1):2 with lines title "Merkin et al (2010)"
