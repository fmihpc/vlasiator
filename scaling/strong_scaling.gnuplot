case = "light"
#case = "medium"
#case = "heavy"

set log x
set xrange [1:513]
set yrange [10:11000]
set log y
set xlabel "Nodes"
set ylabel "Total run time (s)"
set title "Strong scaling (".case.")"
set key right bottom

#set terminal qt
set term png font "Corbel,14"
set output "strong_scaling_".case.".png"


plot "202506_CSC_Mahti/strong/timings_".case.".dat" u 1:2 w lp lw 2 t "Total run time Mahti (s)", \
     t=0 "202506_CSC_Mahti/strong/timings_".case.".dat" u 1:(t==0?y0=$1*$2:y0, t=t+1, y0/$1) w lp lw 2 t "Ideal scaling Mahti", \
     "202507_EuroHPC_LUMI-C/strong/timings_".case.".dat" u 1:2 w lp lw 2 t "Total run time LUMI-C (s)", \
     u=0 "202507_EuroHPC_LUMI-C/strong/timings_".case.".dat" u 1:(u==0?z0=$1*$2:z0, u=u+1, z0/$1) w lp lw 2 t "Ideal scaling LUMI-C"
