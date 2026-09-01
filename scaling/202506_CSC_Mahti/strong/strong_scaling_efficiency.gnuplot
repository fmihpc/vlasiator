#case = "light"
#case = "medium"
case = "heavy"

set log x
set xrange [1:201]
set yrange [0:1.1]
set xlabel "Nodes"
set ylabel "Efficiency"
set title "Strong scaling efficiency ".case." (Mahti)"
set key left bottom

set term png font "Corbel,14"
set output "strong_scaling_efficiency_".case.".png"

plot t=0 "timings_".case.".dat" u 1:(t==0?y0=$2:y0, t==0?nc=$1:nc, t=t+1, y0*nc/($1*$2)) w lp lw 2 t "Total run time", \
     t=0 "timings_".case.".dat" u 1:(t==0?y0=$3:y0, t==0?nc=$1:nc, t=t+1, y0*nc/($1*$3)) w lp lw 2 t "Propagate", \
     t=0 "timings_".case.".dat" u 1:(t==0?y0=$4:y0, t==0?nc=$1:nc, t=t+1, y0*nc/($1*$4)) w lp lw 2 t "Spatial-space", \
     t=0 "timings_".case.".dat" u 1:(t==0?y0=$5:y0, t==0?nc=$1:nc, t=t+1, y0*nc/($1*$5)) w lp lw 2 t "Velocity-space", \
     1 lw 2 lc "black" notitle

