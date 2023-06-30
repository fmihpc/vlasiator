#!/usr/bin/gnuplot
#
# Before running this, do:
#
# make differentialFlux
# for i in 1e5 1e6 1e7 1e8; do ./differentialFlux 1e6 $i > $i.dat; done
#

set terminal pngcairo size 800,600
set output "diffFlux.png"

set title "Differential flux for electron population at n=10^6 m^{-3}"
set xlabel "Energy (keV)"
set ylabel "Differential flux (eV/m^2/s)"

#set logscale x
#set logscale y
#set yrange [1e4:1e15]

#plot '1e5.dat' with histeps title "T = 10^5 K",\
# '1e6.dat' with histeps title "T=10^6 K",\
# '1e7.dat' with histeps title "T=10^7 K",\
# '1e8.dat' with histeps title "T=10^8 K"

set xrange [0:10]
#set yrange [0:1e9]
set yrange [0:]
set grid
set arrow from 1,graph 0 to 1, graph 1 nohead lc rgbcolor "grey"
plot '1kev.dat' with histeps title "T = 1keV"
