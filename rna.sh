#!/bin/bash
#SBATCH --job-name=main
#SBATCH --output=batch_job.log
#SBATCH --time=01:00:00
#SBATCH --account=PAS2967  # Specify PI's account here
#SBATCH --ntasks=1                   # Number of parallel tasks (processes)
#SBATCH --cpus-per-task=48            # CPUs allocated per task (adjust as needed)
#SBATCH --constraint=48core



# Your job's commands go here
source ~/.bashrc
conda activate env1
python rna.py