#case = "light"
#case = "medium"
case = "heavy"

set log x
set xrange [1:201]
set yrange [10:11000]
set log y
set xlabel "Nodes"
set ylabel "Total run time (s)"
set title "Strong scaling ".case." (Mahti)"
set key right bottom

set term png font "Corbel,14"
set output "strong_scaling_".case.".png"

plot "timings_".case.".dat" u 1:2 w lp lw 2 t "Total run time (s)", \
     t=0 "timings_".case.".dat" u 1:(t==0?y0=$1*$2:y0, t=t+1, y0/$1) w lp lw 2 t "Ideal scaling"

