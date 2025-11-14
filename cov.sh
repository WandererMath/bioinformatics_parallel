#!/bin/bash
#SBATCH --job-name=main
#SBATCH --output=batch_job.log
#SBATCH --time=04:00:00
#SBATCH --account=PAS2967  # Specify PI's account here
#SBATCH --ntasks=1                   # Number of parallel tasks (processes)
#SBATCH --cpus-per-task=48            # CPUs allocated per task (adjust as needed)
#SBATCH --constraint=48core



# Your job's commands go here
source ~/.bashrc
conda activate bioinfo
python cov_plot.py