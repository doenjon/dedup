#!/bin/bash

# Run each config in its own job
sbatch run_configs/run_Av.sh
sbatch run_configs/run_Sl.sh
sbatch run_configs/run_Mm.sh
sbatch run_configs/run_Dei.sh
sbatch run_configs/run_At.sh
#sbatch run_configs/run_Pt.sh
