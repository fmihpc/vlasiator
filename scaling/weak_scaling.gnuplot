set log x
set xrange [1:513]
set xlabel "Nodes"
set ylabel "Total run time (s)"
set title "Weak scaling"
set key left top

#set terminal qt
set term png font "Corbel,14"
set output "weak_scaling.png"

plot t=0 "202506_CSC_Mahti/weak/along_y.dat" u 1:(t==0?nc=$2:nc, t=t+1, $2) w lp lw 2 t "Extended along y (Mahti)", \
         "202506_CSC_Mahti/weak/along_z.dat" u 1:2 w lp lw 2 t "Extended along z (Mahti)", \
         "202506_CSC_Mahti/weak/along_yz.dat" u 1:2 w lp lw 2 t "Extended along y, z (Mahti)", \
     u=0 "202507_EuroHPC_LUMI-C/weak/along_y.dat" u 1:(u==0?nc=$2:nc, u=u+1, $2) w lp lw 2 t "Extended along y (LUMI-C)", \
         "202507_EuroHPC_LUMI-C/weak/along_z.dat" u 1:2 w lp lw 2 t "Extended along z (LUMI-C)", \
         "202507_EuroHPC_LUMI-C/weak/along_yz.dat" u 1:2 w lp lw 2 t "Extended along y, z (LUMI-C)", \
         nc lw 2 lc "black" notitle

