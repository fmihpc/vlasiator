#!/bin/bash

# ==============================
# User-configurable parameters
# ==============================

# SLURM job parameters
NODES=$1
NTASKS=$2
GPUS=$3
TASKS_PER_NODE=$(( NTASKS / NODES ))
GPUS_PER_NODE=$(( GPUS / NODES ))

# Grid dimensions
X_LENGTH=$4
Y_LENGTH=$5
Z_LENGTH=$6

X_MIN=$7
X_MAX=$8
Y_MIN=$9
Y_MAX=${10}
Z_MIN=${11}
Z_MAX=${12}

T_MAX=${13}
PARTITION=${14}
FOLDER=${15}
SH_INPUT=${16}

SH_FILE="$FOLDER/weakScaling_${GPUS}/weakScaling_${GPUS}.sh"
CFG="weakScaling_${GPUS}"
CFG_FILE="$FOLDER/weakScaling_${GPUS}/weakScaling_${GPUS}.cfg"

CWD=$(pwd)

mkdir $FOLDER/weakScaling_${GPUS}
cp ./$FOLDER/sw1.dat $FOLDER/weakScaling_${GPUS}

# ==============================
# Generate weak_scaling.sh
# ==============================

cp "$SH_INPUT" "$SH_FILE"

chmod +x ${SH_FILE}

# ==============================
# Generate weakScaling.cfg
# ==============================

awk -v XL="$X_LENGTH" -v YL="$Y_LENGTH" -v ZL="$Z_LENGTH" \
    -v XM="$X_MIN" -v XMx="$X_MAX" -v YM="$Y_MIN" -v YMx="$Y_MAX" \
    -v ZM="$Z_MIN" -v ZMx="$Z_MAX" -v TM="$T_MAX" '
BEGIN {in_section=0}
{
    # detect start/end of [gridbuilder] section
    if ($0 ~ /^\[gridbuilder\]/) {in_section=1; print; next}
    else if ($0 ~ /^\[/) {in_section=0; print; next}

    if (in_section) {
        # replace keys if they match
        if ($1 ~ /^x_length/)  {print "x_length = " XL; next}
        if ($1 ~ /^y_length/)  {print "y_length = " YL; next}
        if ($1 ~ /^z_length/)  {print "z_length = " ZL; next}
        if ($1 ~ /^x_min/)     {print "x_min = " XM; next}
        if ($1 ~ /^x_max/)     {print "x_max = " XMx; next}
        if ($1 ~ /^y_min/)     {print "y_min = " YM; next}
        if ($1 ~ /^y_max/)     {print "y_max = " YMx; next}
        if ($1 ~ /^z_min/)     {print "z_min = " ZM; next}
        if ($1 ~ /^z_max/)     {print "z_max = " ZMx; next}
        if ($1 ~ /^t_max/)     {print "t_max = " TM; next}
    }

    # print other lines unchanged
    print
}' "$FOLDER/$FOLDER.cfg" > "$CFG_FILE"

sed -i "s/__TASKS_PER_NODE__/${TASKS_PER_NODE}/" "$SH_FILE"
sed -i "s/__NTASKS__/${NTASKS}/" "$SH_FILE"
sed -i "s/__GPUS_PER_NODE__/${GPUS_PER_NODE}/" "$SH_FILE"
sed -i "s/__GPUS__/${GPUS}/" "$SH_FILE"
sed -i "s/__NODES__/${NODES}/" "$SH_FILE"
sed -i "s/__PARTITION__/${PARTITION}/" "$SH_FILE"
sed -i "s|__CWD__|${CWD}|g" "$SH_FILE"
sed -i "s/__CFG__/${CFG}/" "$SH_FILE"
sed -i "s/__FOLDER__/${FOLDER}/" "$SH_FILE"

sed -i "s/__X_LENGTH__/${X_LENGTH}/" "$CFG_FILE"
sed -i "s/__Y_LENGTH__/${Y_LENGTH}/" "$CFG_FILE"
sed -i "s/__Z_LENGTH__/${Z_LENGTH}/" "$CFG_FILE"
sed -i "s/__X_MAX__/${X_MAX}/" "$CFG_FILE"
sed -i "s/__Y_MAX__/${Y_MAX}/" "$CFG_FILE"
sed -i "s/__Z_MAX__/${Z_MAX}/" "$CFG_FILE"
sed -i "s/__X_MIN__/${X_MIN}/" "$CFG_FILE"
sed -i "s/__Y_MIN__/${Y_MIN}/" "$CFG_FILE"
sed -i "s/__Z_MIN__/${Z_MIN}/" "$CFG_FILE"
sed -i "s/__T_MAX__/${T_MAX}/" "$CFG_FILE"

echo "✅ Generated weak_scaling.sh and weakScaling/weakScaling.cfg"
