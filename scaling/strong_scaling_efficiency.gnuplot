case = "light"
#case = "medium"
#case = "heavy"

set log x
set xrange [1:513]
set yrange [0:1.1]
set xlabel "Nodes"
set ylabel "Efficiency"
set title "Strong scaling efficiency (".case.")"
set key left bottom

#set terminal qt
set term png font "Corbel,14"
set output "strong_scaling_efficiency_".case.".png"

plot t=0 "202506_CSC_Mahti/strong/timings_".case.".dat" u 1:(t==0?y0=$2:y0, t==0?nc=$1:nc, t=t+1, y0*nc/($1*$2)) w lp lw 2 t "Total run time (Mahti)", \
     t=0 "202506_CSC_Mahti/strong/timings_".case.".dat" u 1:(t==0?y0=$3:y0, t==0?nc=$1:nc, t=t+1, y0*nc/($1*$3)) w lp lw 2 t "Propagate (Mahti)", \
     t=0 "202506_CSC_Mahti/strong/timings_".case.".dat" u 1:(t==0?y0=$4:y0, t==0?nc=$1:nc, t=t+1, y0*nc/($1*$4)) w lp lw 2 t "Spatial-space (Mahti)", \
     t=0 "202506_CSC_Mahti/strong/timings_".case.".dat" u 1:(t==0?y0=$5:y0, t==0?nc=$1:nc, t=t+1, y0*nc/($1*$5)) w lp lw 2 t "Velocity-space (Mahti)", \
     u=0 "202507_EuroHPC_LUMI-C/strong/timings_".case.".dat" u 1:(u==0?z0=$2:z0, u==0?nc=$1:nc, u=u+1, z0*nc/($1*$2)) w lp lw 2 t "Total run time (LUMI-C)", \
     u=0 "202507_EuroHPC_LUMI-C/strong/timings_".case.".dat" u 1:(u==0?z0=$3:z0, u==0?nc=$1:nc, u=u+1, z0*nc/($1*$3)) w lp lw 2 t "Propagate (LUMI-C)", \
     u=0 "202507_EuroHPC_LUMI-C/strong/timings_".case.".dat" u 1:(u==0?z0=$4:z0, u==0?nc=$1:nc, u=u+1, z0*nc/($1*$4)) w lp lw 2 t "Spatial-space (LUMI-C)", \
     u=0 "202507_EuroHPC_LUMI-C/strong/timings_".case.".dat" u 1:(u==0?z0=$5:z0, u==0?nc=$1:nc, u=u+1, z0*nc/($1*$5)) w lp lw 2 t "Velocity-space (LUMI-C)", \
     1 lw 2 lc "black" notitle

