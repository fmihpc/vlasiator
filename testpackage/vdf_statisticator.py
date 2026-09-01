

import numpy as np
import scipy.ndimage as ndimage
import analysator as pt
import matplotlib.pyplot as plt
import os
import scipy.stats

#text file of the run names
#I used run_[timestamp]_[processid], so i made the list with `ls -l | grep -Po "run_\d+_\d+` > runlist.txt
#within those run folders are the test folders
runs = np.loadtxt("/home/siclasse/vlasiator_hile/vlasiator/testpackage/acctest500/runlist.txt", object)

run_title="acctest"
do_slices=True            #plot mean,std,snr vdf in slices not using analysator
do_slices_analysator=False #this should maybe be done with ohters disabled, uses analysator to plot the slices of the vdf instead 
do_integrated=True          #use analysator to get the integrated grid (note that i modified analysator to return the mesh and bins earlier in plot_vdf())
do_peak=True                #get the peak indices of the peak value for the vdf from the first run and make historgram of the vdf value in that cell
do_statistics=True          #histogram of bulkv(vg_v),proton/numrho(proton/vg_rho) and bulkrho (vg_rhom)

i = 0
dit = {}
grid_init = []
all_a = []
max_i = len(runs) #for manually setting how many runs to plot
# max_i = 6
peak_id = -1
peak_vals = []
peak_coord_og = np.array([-1, -1, -1])

#I love eating RAM, nom nom
grid_size=[400,400,400]
if do_integrated:
    bins_runs=np.full((max_i,grid_size[0],grid_size[1]),0,dtype=np.float32)
if do_slices:
    fullgrid = np.full((max_i, *grid_size), 0, dtype=np.float32)

bulkv=[]
bulkrho=np.array([])
numrho=np.array([])

filename="fullf"
for d,run in enumerate(runs):
    ref = f"/home/siclasse/vlasiator_hile/vlasiator/testpackage/acctest500/{run}/acctest_1_maxw_500k_100k_20kms_10deg"
    print(ref)
    if not os.path.isfile(ref+f"/{filename}.0000001.vlsv"):
        print(f"File {ref}/{filename} not found")
        continue
    file = pt.vlsvfile.VlsvReader(ref + f"/{filename}.0000001.vlsv",indexer="dict")
    R_EARTH = 6.371e6
    x = 15
    y = 0
    z = 0

    cellid=1
    # cellid = file.get_cellid_with_vdf([x*R_EARTH, y*R_EARTH, z*R_EARTH])
    if do_statistics:
        bulkv.append(np.array(file.read_variable("vg_v",cellids=cellid)))
        bulkrho=np.append(bulkrho,file.read_variable("vg_rhom",cellids=cellid))
        numrho=np.append(numrho,file.read_variable("proton/vg_rho",cellids=cellid))

    print(run,cellid)
    pop = "proton"
    velocity_cell_map = file.read_velocity_cells(cellid, pop)
    maps = list(zip(*velocity_cell_map.items()))
    velocity_cell_ids = np.array(maps[0], dtype=np.int64)
    velocity_cell_values = np.array(maps[1], dtype=np.float32)

    velocity_cell_coordinates = file.get_velocity_cell_coordinates(
        velocity_cell_ids, pop
    )
    velocity_cell_indices = file.get_velocity_cell_indices(
        velocity_cell_coordinates, pop
    )

    if i == 0:
        dv = (velocity_cell_coordinates[1] - velocity_cell_coordinates[0])[0]
        peak_id = np.argmax(velocity_cell_values)
        peak_coord = velocity_cell_coordinates[peak_id]
        peak_highest_corner = np.max(velocity_cell_indices, axis=0)
        peak_lowest_corner = np.min(velocity_cell_indices, axis=0)
        peak_ind = velocity_cell_indices[peak_id]
        peak_id = velocity_cell_ids[peak_id]
        slicevecs=[[0,0,peak_coord[2]+z*dv] for z in range(-10,10)]
   
    #note that the analysator plot_vdf function was modified to return these values after they are set the final time
    if do_integrated:
        xmesh,ymesh,bins=pt.plot.plot_vdf(vlsvobj=file,cellids=[cellid],outputfile=f'./plot_vdf_{d}.png',xy=1,slicethick=1,fmin=1e-15,fmax=4e-9,box=np.array([-500,500,-1000,1000])*10**3)
        bins_runs[i]=bins
        plt.clf()

    if do_slices_analysator:
        for slicevec in slicevecs:
            os.system(f"mkdir slices_{slicevec[2]}_analysator_plot")
            pt.plot.plot_vdf(vlsvobj=file,cellids=[1],outputfile=f'./slices_{slicevec[2]}_analysator_plot/plot_vdf_{slicevec[2]}_{d}_{run_title}.png',center=slicevec,xy=1,slicethick=0.5,box=np.array([-500,500,-1000,1000])*10**3,reducer="average",fmin=1e-15,fmax=4e-9)
    
    if do_peak:
        peak_vals.append(
            velocity_cell_values[np.argwhere(velocity_cell_ids == peak_id)][0][0]
        )

    if do_slices:
        for id, cellindices in enumerate(velocity_cell_indices):
            ind = np.asarray(cellindices)
            fullgrid[i, ind[0],ind[1],ind[2]] = velocity_cell_values[id]

    i += 1
    if i == max_i:
        break



if do_statistics:
    bulkv=np.array(bulkv)
    for name,val in [("numrho",numrho),("bulkv_x",bulkv[:,0]),("bulkv_y",bulkv[:,1]),("bulkv_z",bulkv[:,2]),("bulkrho",bulkrho*10**25)]:
        plt.title(f"{run_title} {name} hist")

        statistic,pval=scipy.stats.shapiro(val)
        plt.hist(val, bins="auto", density=False,label=f"shapiro-wilk={statistic:.3f},pval={pval:.3f},mean={np.mean(val):.3f},std={np.std(val,ddof=1):.3f}")
        #if chosen alpha level (the level at which we are comfortable falsely rejecting null(data is normally distributed) is above p-value, we may not reject the possibility that it is normally distributed) 
        #note that this is not the same as to say it is normally distributed or not.

        #had to scale with e-25 for some so the shapiro-wilk didn't complain SHOULD BE CHANGED IF USING DIFFERENT DATA
        
        plt.gca().set_xlabel("bins e-25" if name=="bulkrho" else "bins")
        plt.gca().set_ylabel("counts")
        plt.gcf().set_size_inches(8,6)
        plt.legend()
        plt.savefig(f"{name}_hist_{run_title}.png")
        plt.clf()

if do_peak:
    plt.hist(peak_vals, bins="auto", density=False)
    plt.savefig(f"./peak_hist_{run_title}_{i}_peakid_{peak_id}.png")
    plt.clf()


if do_slices:
    std = np.nanstd(fullgrid, axis=0, ddof=1) #should not have nans but whatever
# std = np.nan_to_num(std)

    mean = np.nanmean(fullgrid, axis=0)
    print("std maximum ",np.max(std))
    tick_n = 10
    xtick = (
        np.arange(0, grid_size[0], 1) * dv  + file.get_velocity_mesh_extent()[0:3][0] + dv 
    )
    xtick = [f"{tick * 10**-3:.0f}" for tick in xtick]
    snr = np.full(std.shape, np.nan)
    np.divide(mean, std, where=std > 0, out=snr)

# snr_peak_coords=[]
# peak_snr=np.array([-1,-1,-1])

if do_integrated:
    mean_inte=np.mean(bins_runs,axis=0)
    std_inte=np.std(bins_runs,axis=0)

    snr_inte = np.full(std_inte.shape,0,dtype=np.float32)
    np.divide(mean_inte, std_inte, where=std_inte > 0, out=snr_inte)

for scaling in ["linear", "log"]:
        if do_integrated: 
            for name, data in [("std", std_inte), ("mean", mean_inte), ("snr", snr_inte)]:
                plt.gcf().set_size_inches(8,8)
                plt.gca().set_xlim(-800,800)
                plt.gca().set_ylim(-800,800)
                plt.pcolormesh(xmesh,ymesh,data,norm=scaling)

                plt.colorbar(orientation="horizontal")
                plt.title(
                    f"{name} of VDF over {i} runs, scaling {scaling}, run {run_title}"
                )

                fig = plt.gcf()
                ax = plt.gca()

                plt.xlabel("v_x [m/s]")
                plt.xticks(rotation=45)
                plt.ylabel("v_y [m/s]")
                plt.savefig(f"./integrated_{name}_{i}_{scaling}_{run_title}.png", dpi=300)
                ax.clear()
                fig.clear()   

        if do_slices:
            for slice_i in range(peak_ind[2]-10, peak_ind[2]+10):
                for name, data in [("std", std), ("mean", mean), ("snr", snr)]:
                    if np.max(data[:, ::-1, slice_i].T) != np.min(data[:, ::-1, slice_i].T):
                        cutoff = 1 * (mean[:, :, slice_i] > 10 ** (-16))

                        x, y = np.meshgrid(
                            np.arange(0, grid_size[0]), np.arange(0, grid_size[1])
                        )
                        plt.contour(
                            x,
                            y,
                            cutoff.T,
                            0,
                            colors="black",
                            linewidths=0.5,
                            linestyles="dashed",
                        )

                        #following was for tracking the SNR blob
                        # if name=="snr" and scaling=="linear":
                        #     peak_snr=np.argwhere(data[:,:,slice_i].T==np.nanmax(data[:,:,slice_i].T))[0]
                        #     snr_peak_coords.append([peak_snr[0],peak_snr[1],slice_i])
                        #     continue
                        # else:
                        #     continue
                        #

                        plt.imshow(data[:, :, slice_i].T, norm=scaling)

                        plt.colorbar(orientation="horizontal")
                        plt.title(
                            f"{name} of VDF over {i} runs, scaling {scaling} slice z={slice_i}, run {run_title}"
                        )

                        fig = plt.gcf()
                        ax = plt.gca()
                        ax.set_xticks(
                            np.linspace(0.5, grid_size[0] - 0.5, grid_size[0])[::tick_n]
                        )
                        ax.set_yticks(
                            np.linspace(-0.5, grid_size[1] - 1.5, grid_size[1])[::-tick_n]
                        )
                        ax.set_xticklabels(xtick[::tick_n])
                        ax.set_yticklabels(xtick[::-tick_n])
                        fig.set_size_inches(6, 6)

                        fig.tight_layout()

                        plt.xlim(peak_lowest_corner[0] - 10, peak_highest_corner[0] + 10)
                        plt.ylim(peak_lowest_corner[1], peak_highest_corner[1])
                        plt.xlabel("v_x [km/s]")
                        plt.xticks(rotation=45)
                        plt.ylabel("v_y [km/s]")
                        plt.savefig(f"{name}_{slice_i}_{i}_{scaling}_{run_title}.png", dpi=300)
                        ax.clear()
                        fig.clear()

