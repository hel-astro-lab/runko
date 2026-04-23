#!/bin/bash

if [ $SLURM_PROCID -eq 1 ]; then
	rocprof-sys-python-3.11 -- ../pic.py
else 
	python ../pic.py
fi
