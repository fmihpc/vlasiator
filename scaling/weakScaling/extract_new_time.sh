#!/bin/bash
# Usage: ./extract_new_time.sh file.out

FILE=$1
PERCENT_DIFF=$2

# Find the line with "NA" and extract the 2nd column (new_time)
new_time=$(awk -F'|' '/NA[ \t]*\|/ {gsub(/[ \t]/,"",$2); print $2}' "$FILE")

# Compute scaled time assuming percent diff = 0
scaled_time=$(awk -v t="$new_time" -v p="$PERCENT_DIFF" 'BEGIN {
    scaled = t / (1 + p/100.0)
    print scaled
}')

echo "$scaled_time" >> results.txt
