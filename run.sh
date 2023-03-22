#!/bin/bash

#SBATCH --job-name="TDSE3.9"
#SBATCH --export=ALL
#SBATCH --nodes=1
#SBATCH --cpus-per-task=96
srun python3 launch.py
