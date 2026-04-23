#!/bin/bash

if [ $SLURM_PROCID -eq 1 ]; then
	rocprofv3 --stats --hip-trace --output-format pftrace -- python ../pic.py
else 
	python ../pic.py
fi
