#!/bin/bash

if [ "$#" -ne 2 ]; then
   echo "usage: $0 path-to-logfile run_id"
   exit 1
fi

LOGF=$1
RUN_ID=$2

   # Vlasiator iowrite.cpp outputs for reference
   # logFile << "(writeGrid) Wrote ";

   # if (bytesWritten > 1.0e9) logFile << bytesWritten/1.0e9 << " GB in ";
   # else if (bytesWritten > 1e6) logFile << bytesWritten/1.0e6 << " MB in ";
   # else if (bytesWritten > 1e3) logFile << bytesWritten/1.0e3 << " kB in ";
   # else logFile << bytesWritten << " B in ";

   # logFile << writeTime << " seconds, approximate data rate is ";

   # if (bytesWritten/writeTime > 1e9) logFile << bytesWritten/writeTime/1e9 << " GB/s";
   # else if (bytesWritten/writeTime > 1e6) logFile << bytesWritten/writeTime/1e6 << " MB/s";
   # else if (bytesWritten/writeTime > 1e3) logFile << bytesWritten/writeTime/1e3 << " kB/s";
   # else logFile << bytesWritten/writeTime << " B/s";
   # logFile << endl;


grep '(writeGrid) Wrote' $LOGF > "$RUN_ID"_writerates.txt

awk 'BEGIN {print("[GB]\t[s]\t[GB/s]")} {if($4=="GB"){ print($3 "\t" $6 "\t" $12) } else if($4=="MB"){writesize=$3/1000; rate=$12/1000; print(writesize "\t" $6 "\t" rate) } else if($3=="kB"){writesize=$3/10000000; rate=$12/10000000; print(writesize "\t" $6 "\t" rate) }}' "$RUN_ID"_writerates.txt > "$RUN_ID"_rates.txt

python plot_rates.py $RUN_ID