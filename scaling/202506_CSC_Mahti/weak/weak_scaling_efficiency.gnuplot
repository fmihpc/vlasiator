set log x
set xrange [1:201]
set yrange [0:1.1]
set xlabel "Nodes"
set ylabel "Efficiency"
set title "Weak scaling efficiency (Mahti)"
set key left bottom

set term png font "Corbel,14"
set output "weak_scaling_efficiency.png"

plot t=0 "all.dat" u 1:(t==0?y0=$2:y0, t=t+1,y0/$2) w p lw 2 t "Total run time" , \
     t=0 "all.dat" u 1:(t==0?y0=$3:y0, t=t+1,y0/$3) w p lw 2 t "Propagate", \
     t=0 "all.dat" u 1:(t==0?y0=$4:y0, t=t+1,y0/$4) w p lw 2 t "Spatial-space", \
     t=0 "all.dat" u 1:(t==0?y0=$5:y0, t=t+1,y0/$5) w p lw 2 t "Velocity-space", \
     1 w l lw 2 lc "black" notitle



