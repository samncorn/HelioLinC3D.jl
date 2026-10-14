#!/bin/bash

#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --job-name=heliolinc-optimization
#SBATCH --account=samuel41-ic

#SBATCH --output=jobarray.o%A_%a
#SBATCH --error=jobarray.e%A_%a

#SBATCH --mail-user=samuel41@illinois.edu
#SBATCH --mail-type=BEGIN,END

export HELIOLINC_DIR=""

# source clustering parameter from job array

cd # project dir

detections_dir=
observer_file=
hypothesis_file=

night=15

# couldn't get the cli wrapper to work on the cluster so we have to manually specify the julia program to run
julia -t 128 "$HELIOLINC_DIR/scripts/heliolinc_v2.jl" \
    --detections="$detections_dir" \
    --observer="$observer_file" \
    --hypothesis="$hypothesis_file" \
    --night "$night" \
    --out "clusters/$night/" \
    --p 128
