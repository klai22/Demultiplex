#!/bin/bash

#sbatch demultiplex.sh

#SBATCH --job-name=demultiplexp3
#SBATCH --account=bgmp
#SBATCH --partition=gpu
#SBATCH --cpus-per-task=5 
#SBATCH --mem=100G
#SBATCH --nodes=1
#SBATCH --output=/home/kenlai/bgmp/bioinfo/Bi622/Assignments/Demultiplex/Assignment-the-third/o_e_files/slurm-special-name-%j.out
#SBATCH --error=/home/kenlai/bgmp/bioinfo/Bi622/Assignments/Demultiplex/Assignment-the-third/o_e_files/slurm-special-name-%j.err

/usr/bin/time -v ./demultiplex.py -R2qs 29 -R3qs 26 -hd 3 -R1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz -R2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz -R3 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz -R4 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz -bars /projects/bgmp/shared/2017_sequencing/indexes.txt

exit