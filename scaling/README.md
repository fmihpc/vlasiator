# Scaling results obtained on a variety of architectures

## Architectures

### CSC's Mahti
[CSC's Mahti](https://docs.csc.fi/computing/systems-mahti/) CPU partition, 2x AMD Rome 7H12 CPUs per node at 2.6 GHz base frequency, 3.3 GHz max boost. 16 tasks x 16 threads per node (i.e. multithreading).

### EuroHPC's LUMI-C
[EuroHPC's LUMI-C](https://docs.lumi-supercomputer.eu/), 2x AMD EPYC 7763 CPUs per node at 2.45 GHz base frequency, 3.5 GHz max boost. 16 tasks x 16 threads per node (i.e. multithreading).


## Setup
Identical IPShock-based setups were run in summer 2025 on [CSC's Mahti](https://docs.csc.fi/computing/systems-mahti/) and [EuroHPC's LUMI-C](https://docs.lumi-supercomputer.eu/), see the subfolders. Those were run on dev with slightly updated version for the latter.

[CSC/Mahti results](202506_CSC_Mahti/README.md)

[EuroHPC/LUMI-C results](202507_EuroHPC_LUMI-C/README.md)


## Plots
Plot the results with `gnuplot -p <script>.gnuplot` for png plots, comment the relevant lines to obtain interactive plots instead.

### Weak scaling
Combined plots of weak scaling and weak scaling efficiency from [CSC/Mahti](202506_CSC_Mahti/weak/) and [EuroHPC/LUMI-C](202507_EuroHPC_LUMI-C/weak/).

![](weak_scaling.png)
![](weak_scaling_efficiency.png)

### Strong scaling

Combined plots of strong scaling and strong scaling efficiency from [CSC/Mahti](202506_CSC_Mahti/weak/) and [EuroHPC/LUMI-C](202507_EuroHPC_LUMI-C/weak/).

Three cases run, "light", "medium", and "heavy" with progressively larger problems (extended in *y* and *z* and lower phase-space density threshold).

![](strong_scaling_light.png)
![](strong_scaling_medium.png)
![](strong_scaling_heavy.png)

![](strong_scaling_efficiency_light.png)
![](strong_scaling_efficiency_medium.png)
![](strong_scaling_efficiency_heavy.png)
