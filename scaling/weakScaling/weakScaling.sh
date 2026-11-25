#!/bin/bash
# ==============================
# Launcher for weak scaling
# ==============================

# User-defined parameters

SCALE_X=1
SCALE_Y=0
SCALE_Z=0

PARTITION="standard-g"

### USER EDITED SECTION ENDS

# Input arguments

FOLDER=$1
SH_INPUT=$2
MIN_EXPONENT=$3
MAX_EXPONENT=$4

### INPUT ARGUMENTS END

DIMENSION=$(( SCALE_X + SCALE_Y + SCALE_Z ))

# Initialize variables
X_LENGTH_BASE=0
Y_LENGTH_BASE=0
Z_LENGTH_BASE=0
X_MIN_BASE=0
X_MIN_EXP=0
X_MAX_BASE=0
X_MAX_EXP=0
Y_MIN_BASE=0
Y_MIN_EXP=0
Y_MAX_BASE=0
Y_MAX_EXP=0
Z_MIN_BASE=0
Z_MIN_EXP=0
Z_MAX_BASE=0
Z_MAX_BASE=0
T_MAX_BASE=0
T_MAX_EXP=0

# Function to split scientific notation and convert to float
split_sci() {
    val="$1"
    if [[ "$val" =~ ^(-?[0-9.]+)[eE]([+-]?[0-9]+)$ ]]; then
        base="${BASH_REMATCH[1]}"
        exp="${BASH_REMATCH[2]}"
    else
        base="$val"
        exp=0
    fi
    echo "$base $exp"
}

# Flag for section
in_section=0

while IFS= read -r line || [[ -n "$line" ]]; do
    # Trim spaces
    line="${line#"${line%%[![:space:]]*}"}"
    line="${line%"${line##*[![:space:]]}"}"
    [[ -z "$line" ]] && continue
    [[ "$line" =~ ^# ]] && continue

    # Start/end of section
    if [[ "$line" =~ ^\[[[:space:]]*gridbuilder[[:space:]]*\] ]]; then
        in_section=1
        continue
    elif [[ "$line" =~ ^\[.*\] ]]; then
        in_section=0
        continue
    fi

    # Process lines inside section
    if [[ $in_section -eq 1 ]] && [[ "$line" == *=* ]]; then
        key=${line%%=*}
        val=${line#*=}
        key=$(echo "$key" | tr -d ' ')
        val=$(echo "$val" | tr -d ' ')
        
        # Optional: get base and exponent
        read base exp <<< $(split_sci "$val")
        
        # Update corresponding Bash variable
        case "$key" in
            x_length) X_LENGTH_BASE=$base ;;
            y_length) Y_LENGTH_BASE=$base ;;
            z_length) Z_LENGTH_BASE=$base ;;
            x_min)    X_MIN_BASE="$base";    X_MIN_EXP="$exp" ;;
            x_max)    X_MAX_BASE="$base";    X_MAX_EXP="$exp" ;;
            y_min)    Y_MIN_BASE="$base";    Y_MIN_EXP="$exp" ;;
            y_max)    Y_MAX_BASE="$base";    Y_MAX_EXP="$exp" ;;
            z_min)    Z_MIN_BASE="$base";    Z_MIN_EXP="$exp" ;;
            z_max)    Z_MAX_BASE="$base";    Z_MAX_EXP="$exp" ;;
            t_max)    T_MAX_BASE="$base";    T_MAX_EXP="$exp" ;;
            # Add other keys if needed
        esac
    fi
done < "$FOLDER/$FOLDER.cfg"


: > percent_diff.txt
# ==============================
# Execute the generator script
# ==============================
# Assume the generator script is called `generate_weak_scaling.sh`
# and is in the same folder
for (( i=MIN_EXPONENT; i<=MAX_EXPONENT; i++ )); do
    j=$(( i - MIN_EXPONENT ))
    NTASKS=$(echo "2^$i" | bc -l)
    GPUS=$NTASKS
    NODES=$(( ($GPUS + 8 - 1) / 8 ))

    X_LENGTH=$(echo "$j $X_LENGTH_BASE $SCALE_X $DIMENSION" | awk '{printf "%d\n", ((2^($3/$4))^$1)*$2}')
    Y_LENGTH=$(echo "$j $Y_LENGTH_BASE $SCALE_Y $DIMENSION" | awk '{printf "%d\n", ((2^($3/$4))^$1)*$2}')
    Z_LENGTH=$(echo "$j $Z_LENGTH_BASE $SCALE_Z $DIMENSION" | awk '{printf "%d\n", ((2^($3/$4))^$1)*$2}')

    X_MIN=$(echo "($X_LENGTH * $X_MIN_BASE)/$X_LENGTH_BASE" | bc -l)e${X_MIN_EXP}
    Y_MIN=$(echo "($Y_LENGTH * $Y_MIN_BASE)/$Y_LENGTH_BASE" | bc -l)e${Y_MIN_EXP}
    Z_MIN=$(echo "($Z_LENGTH * $Z_MIN_BASE)/$Z_LENGTH_BASE" | bc -l)e${Z_MIN_EXP}
    X_MAX=$(echo "($X_LENGTH * $X_MAX_BASE)/$X_LENGTH_BASE" | bc -l)e${X_MAX_EXP}
    Y_MAX=$(echo "($Y_LENGTH * $Y_MAX_BASE)/$Y_LENGTH_BASE" | bc -l)e${Y_MAX_EXP}
    Z_MAX=$(echo "($Z_LENGTH * $Z_MAX_BASE)/$Z_LENGTH_BASE" | bc -l)e${Z_MAX_EXP}
    T_MAX=${T_MAX_BASE}e${T_MAX_EXP}

    percent_diff=$(awk -v XL="$X_LENGTH" -v YL="$Y_LENGTH" -v ZL="$Z_LENGTH" \
                        -v XLB="$X_LENGTH_BASE" -v YLB="$Y_LENGTH_BASE" -v ZLB="$Z_LENGTH_BASE" \
                        -v i="$j" 'BEGIN {
        current = XL * YL * ZL
        base_scaled = (2^i) * XLB * YLB * ZLB
        diff = (current - base_scaled) / base_scaled * 100
        print diff
    }')
    echo "$GPUS $percent_diff" >> "percent_diff.txt"

    ./generate_weak_scaling.sh "$NODES" "$NTASKS" "$GPUS" \
                               "$X_LENGTH" "$Y_LENGTH" "$Z_LENGTH" \
                               "$X_MIN" "$X_MAX" "$Y_MIN" "$Y_MAX" \
                               "$Z_MIN" "$Z_MAX" "$T_MAX" "$PARTITION" \
                               "$FOLDER" "$SH_INPUT"
    SH_FILE="$FOLDER/weakScaling_${GPUS}/weakScaling_${GPUS}.sh"
    sbatch "$SH_FILE"

    while true; do
       # Count how many jobs the user has in the given partition
       job_count=$(squeue -u $USER -p dev-g | tail -n +2 | wc -l)

        if [ "$job_count" -lt 2 ]; then
            # Less than 2 jobs, exit the loop
            break
        fi

        # Otherwise, print current jobs and wait
        echo "You currently have $job_count jobs in partition $PARTITION. Waiting..."
        squeue -u $USER -p dev-g
        sleep 30  # wait 10 seconds before checking again
    done

    sleep 5
done
