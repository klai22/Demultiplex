#!/bin/bash

#sbatch demux_stats.sh

#SBATCH --job-name=demux_stats
#SBATCH --account=bgmp
#SBATCH --partition=gpu
#SBATCH --cpus-per-task=5 
#SBATCH --mem=100G
#SBATCH --nodes=3
#SBATCH --output=/home/kenlai/bgmp/bioinfo/Bi622/Assignments/Demultiplex/Assignment-the-third/o_e_files/slurm-special-name-%j.out
#SBATCH --error=/home/kenlai/bgmp/bioinfo/Bi622/Assignments/Demultiplex/Assignment-the-third/o_e_files/slurm-special-name-%j.err

#activate conda env w/ matplotlib 
conda activate bgmp_py312

#run script 

/usr/bin/time -v ./demux_stats.py

exit