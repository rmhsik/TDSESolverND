#!/bin/bash

#SBATCH --job-name="TDSE"
#SBATCH --export=ALL
#SBATCH --nodes=1
#SBATCH --nodelist=nodo08
#SBATCH --cpus-per-task=90
srun venv/bin/python3 launch.py
