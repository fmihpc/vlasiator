#!/usr/bin/gnuplot   
#
# Before running this, do:
# 
# make sigmaProfiles
# ./sigmaProfiles 1e6 1e7 > atmosphere.dat
#

set terminal pdfcairo size 6,5
set output "atmosphere.pdf"

set multiplot

set logscale x
set xlabel "Conductivity (S/m)"
set ylabel "Altitude (km)"
set format x "10^{%T}"

set size 1,0.5

plot [1e-10:1e5] [0:300] 'atmosphere.dat' using 3:1 with lines title "Pedersen conductivity", '' using 4:1 with lines title "Hall conductivity", '' using 5:1 with lines title "Parallel conductivity"

set origin 0,.5
set xlabel "Electron density (m^{-3})"

plot [1e9:1e13] [0:300] 'atmosphere.dat' using 2:1 with lines title "Free electron density"
