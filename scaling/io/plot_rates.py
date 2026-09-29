import numpy as np
import matplotlib.pyplot as plt
from labellines import *
import sys

if len(sys.argv) != 2:
   print("Usage: python plot_write_rates.py run_id")
else:
   run_id = sys.argv[1]

fn = run_id+"_auto_rates.txt"

data = np.loadtxt(fn, skiprows=1)

# Filter out anything below 1s duration or 1GB writesize
# These are assumed failed write, of which there were some
data = data[data[:,1]>1,:]
data = data[data[:,0]>1,:]

peak = np.max(data[:,0]/data[:,1])
print("Maximum performance", peak)
# Reference:
# LUMI-F aggregate bandwidth is 1740 GB/s
# LUMI-P aggregate bandwidth is 240 GB/s

fig, ax = plt.subplots(1)

mapp=ax.scatter(data[:,1],data[:,0], c=data[:,0]/data[:,1], s = 3, cmap='turbo', norm='log')
plt.colorbar(mapp, label="Write performance (GB/s)")
ax.set_xscale("log")
ax.set_yscale("log")
ax.grid()
ax.set_xlabel("time/s")
ax.set_ylabel("Filesize/GB")
fig.draw_without_rendering()

ax.autoscale(False, axis="both")

linelabels = []
for datarate in np.logspace(-3,5,num=17):
   # print(datarate, "GB/s")
   ts = np.linspace(0,2000)
   fs = datarate*ts
   label = "{:.1f} GB/s".format(datarate)
   plt.plot(ts,fs, linestyle=":",color='grey',label=label,linewidth=1)

labelLines(ax.get_lines(),xvals=10,size=6)

ax.set_title('Vlasiator write performance "'+run_id+'"')

plt.savefig(run_id+"_rates.png", dpi=300)
plt.close()

