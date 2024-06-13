#!/usr/bin/gnuplot
#
# Before running this, do:
#
# make sigmaProfiles
# ./sigmaProfiles 1e6 1e7 > atmosphere.dat
#

#set terminal pdfcairo size 4.8,8 dashed
#set output "atmosphere.pdf"
set terminal pngcairo size 960,1600 dashed
set output "atmosphere.png"

set multiplot

set logscale x
set xlabel "Conductivity (S/m)"
set ylabel "Altitude (km)"
set format x "10^{%T}"

set size 1,0.33

plot [1e-10:1e5] [0:250] 'atmosphere.dat' using 3:1 with lines lc 1 dt 1 title "Pedersen conductivity",\
                         '' using 4:1 with lines lc 1 dt 2 title "Hall conductivity",\
                         '' using 5:1 with lines lc 1 dt 3 title "Parallel conductivity",\
                         'atmosphere10.dat' using 3:1 with lines lc 2 dt 1 title "T=10keV",\
                         '' using 4:1 with lines lc 2 dt 2 notitle,\
                         '' using 5:1 with lines lc 2 dt 3 notitle,\
                         'atmosphere100.dat' using 3:1 with lines lc 3 dt 1 title "T=100keV",\
                         '' using 4:1 with lines lc 3 dt 2 notitle,\
                         '' using 5:1 with lines lc 3 dt 3 notitle

set origin 0,.33
set xlabel "Electron density (m^{-3})"

plot [1e9:1e13] [0:250] 'atmosphere.dat' using 2:1 with lines title "Free electron density, T=1kev",\
                        'atmosphere10.dat' using 2:1 with lines title "T=10keV",\
                        'atmosphere100.dat' using 2:1 with lines title "T=100keV"

set origin 0,.66
set xlabel "Production rate (m^{-3} s^{-1})"
set title "Ionosphere with incoming maxwellian electrons n=10^6 m^{-3}"
plot [1e8:1e13] [0:250] 'atmosphere.dat' using 6:1 with lines title "Electron production rate. T=1keV",\
                        'atmosphere10.dat' using 6:1 with lines title "T=10keV",\
                        'atmosphere100.dat' using 6:1 with lines title "T=100keV"
