#!/usr/bin/gnuplot   
#
# Before running this, do:
# 
# make sigmaProfiles
# ./sigmaProfiles 1e6 1e7 > atmosphere.dat
#

set title "Ionosphere with incoming maxwellian electrons n=10^6 m^{-3}, T=1keV"
set terminal pdfcairo size 6,10
set output "atmosphere.pdf"

set multiplot

set logscale x
set xlabel "Conductivity (S/m)"
set ylabel "Altitude (km)"
set format x "10^{%T}"

set size 1,0.33

plot [1e-10:1e5] [0:300] 'atmosphere.dat' using 3:1 with lines title "Pedersen conductivity", '' using 4:1 with lines title "Hall conductivity", '' using 5:1 with lines title "Parallel conductivity"

set origin 0,.33
set xlabel "Electron density (m^{-3})"

plot [1e9:1e13] [0:300] 'atmosphere.dat' using 2:1 with lines title "Free electron density"

set origin 0,.66
set xlabel "Production rate (m^{-3} s^{-1})"
plot [1e8:1e13] [0:300] 'atmosphere.dat' using 6:1 with lines title "Electron production rate"
