set log x
set xrange [1:513]
#set yrange [0:1.1]
set xlabel "Nodes"
set ylabel "Efficiency"
set title "Weak scaling efficiency"
set key left bottom

#set terminal qt
set term png font "Corbel,14"
set output "weak_scaling_efficiency.png"

plot t=0 "202506_CSC_Mahti/weak/all.dat" u 1:(t==0?y0=$2:y0, t=t+1,y0/$2) w p lw 2 ps 3 pt 1 t "Total run time (Mahti)" , \
     t=0 "202506_CSC_Mahti/weak/all.dat" u 1:(t==0?y0=$3:y0, t=t+1,y0/$3) w p lw 2 ps 3 pt 2 t "Propagate (Mahti)", \
     t=0 "202506_CSC_Mahti/weak/all.dat" u 1:(t==0?y0=$4:y0, t=t+1,y0/$4) w p lw 2 ps 3 pt 3 t "Spatial-space (Mahti)", \
     t=0 "202506_CSC_Mahti/weak/all.dat" u 1:(t==0?y0=$5:y0, t=t+1,y0/$5) w p lw 2 ps 3 pt 4 t "Velocity-space (Mahti)", \
     u=0 "202507_EuroHPC_LUMI-C/weak/all.dat" u 1:(u==0?z0=$2:z0, u=u+1,z0/$2) w p lw 2 ps 3 pt 6 t "Total run time (LUMI-C)" , \
     u=0 "202507_EuroHPC_LUMI-C/weak/all.dat" u 1:(u==0?z0=$3:z0, u=u+1,z0/$3) w p lw 2 ps 3 pt 8 t "Propagate (LUMI-C)", \
     u=0 "202507_EuroHPC_LUMI-C/weak/all.dat" u 1:(u==0?z0=$4:z0, u=u+1,z0/$4) w p lw 2 ps 3 pt 10 t "Spatial-space (LUMI-C)", \
     u=0 "202507_EuroHPC_LUMI-C/weak/all.dat" u 1:(u==0?z0=$5:z0, u=u+1,z0/$5) w p lw 2 ps 3 pt 12 t "Velocity-space (LUMI-C)", \
     1 w l lw 2 lc "black" notitle



