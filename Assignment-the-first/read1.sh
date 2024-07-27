#!/bin/bash

#sbatch read1.sh

#SBATCH --job-name=demultiplexp1
#SBATCH --account=bgmp
#SBATCH --partition=gpu
#SBATCH --cpus-per-task=5 
#SBATCH --mem=100G
#SBATCH --nodes=8
#SBATCH --output=/home/kenlai/bgmp/bioinfo/Bi622/Assignments/Demultiplex/Assignment-the-first/o_e_files/slurm-special-name-%j.out
#SBATCH --error=/home/kenlai/bgmp/bioinfo/Bi622/Assignments/Demultiplex/Assignment-the-first/o_e_files/slurm-special-name-%j.err


#Activate matplotlib 
conda activate bgmp_py312

#Running part1.py
/usr/bin/time -v ./part1.py -rl 101 -f /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz -l read1
# /usr/bin/time -v ./part1.py -rl 101 -f /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz -l read2
# /usr/bin/time -v ./part1.py -rl 8 -f /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz -l index1
# /usr/bin/time -v ./part1.py -rl 8 -f /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz -l index2


exit
