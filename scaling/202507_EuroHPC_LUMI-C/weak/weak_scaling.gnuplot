set log x
set xrange [1:513]
set xlabel "Nodes"
set ylabel "Total run time (s)"
set title "Weak scaling (LUMI-C)"
set key left top

set term png font "Corbel,14"
set output "weak_scaling.png"

plot t=0 "along_y.dat" u 1:(t==0?nc=$2:nc, t=t+1, $2) w p lw 2 t "Extended along y", \
         "along_z.dat" u 1:2 w p lw 2 t "Extended along z", \
         "along_yz.dat" u 1:2 w p lw 2 t "Extended along y, z", \
         nc lw 2 lc "black" notitle

