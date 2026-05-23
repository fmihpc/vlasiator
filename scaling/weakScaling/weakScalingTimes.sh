#!/bin/bash

line_index=-1
warn_flag=0  # 0 = no warning, 1 = warning triggered
: > results.txt

while read -r gpus perc_diff; do
    # Pass the percentage difference as an argument
    ./extract_new_time.sh "weak_scaling_${gpus}.out" "$perc_diff"

    # Check if perc_diff is nonzero (allow for floats)
    if [[ $(echo "$perc_diff != 0" | bc -l) -eq 1 ]]; then
        warn_flag=1
    fi
done < ./percent_diff.txt

cat results.txt

# Print a warning if needed
if [[ $warn_flag -eq 1 ]]; then
    echo "WARNING: Using approximate weak scaling with adjusted times!" >&2
fi

