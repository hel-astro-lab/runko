#!/bin/bash

DONE_FLAG="prof_done_${SLURM_JOB_ID}"

if [ $SLURM_PROCID -eq 0 ]; then
    rm -f "$DONE_FLAG"
    rocprof-compute profile -n runko_omniperf --device 0 -- $(which python) ../pic.py
    touch "$DONE_FLAG"
else
    while [ ! -f "$DONE_FLAG" ];
    do
        echo "Rank $SLURM_PROCID starting pass..."
        python ../pic.py
        
        if [ -f "$DONE_FLAG" ]; then
            break
        fi
    done
    echo "Rank $SLURM_PROCID detected completion. Exiting."
fi


sleep 30

if [ $SLURM_PROCID -eq 1 ]; then
    rm -f "$DONE_FLAG"
fi
