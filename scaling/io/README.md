# Helper functions to extract and plot IO performance from Vlasiator logfiles

## get_write_perfs.sh

This script ingests a logfile (and a run ID specifier for temp files and plot outputs) and extracts `(writeGrid)` outputs, converted to GB or GB/s. It then uses `plot_rates.py` to 

## plot_rates.py

This file plots the `get_write_perfs.sh` outputs of write perfomance on a time-filesize scatterplot. Includes a crude filter against failed writes and cannot scrape node counts.

`pip install -r requirements` to get dependencies (`numpy`, `matplotlib`, `matplotlib-label-lines`, of which the last for the snazzy bandwidth labels).

