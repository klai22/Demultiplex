#!/usr/bin/env python

#BFORE RUNNING SCRIPT BE SURE TO ACTIVATE CONDA ENV w/ MATPLOTLIB: conda activate bgmp_py312

#./part1.py -rl 101 -f R1_test.fq.gz -o test.png
#./part1.py -rl 101 -f R1_test.fq.gz -l test

#./part1.py -rl 101 -f /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz -l read1
#./part1.py -rl 101 -f /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz -l read2
#./part1.py -rl 8 -f /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz -l index1
#./part1.py -rl 8 -f /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz -l index2


#setting get_args 
import argparse
#defining input args 
def get_args():
    parser = argparse.ArgumentParser(description="setting global vars for parsing")
    parser.add_argument("-f")
    parser.add_argument("-rl")
    #parser.add_argument("-o")
    parser.add_argument("-l")
    return parser.parse_args()

#setting input args--> variables 
args=get_args()
filename=args.f
read_len=int(args.rl)
#output_file=args.o
label=args.l
# organism=args.org
# biom_file=args.biom
# Afiltered_output_file=args.o2
# Bfiltered_output_file=args.o3
# pid_len=int(args.pl)


#Initialize Empty List
def init_list(lst: list, value: float=0.0) -> list:
    '''This function takes an empty list and will populate it with
    the value passed in "value". If no value is passed, initializes list
    with 101 values of 0.0.'''
    for i in range(read_len):
        lst.append(value)
    return lst


#Importing modules
import bioinfo
import gzip
import matplotlib.pyplot as plt


# Populating List w/ AVG QSs 

#POPULATING w/ QS SUMs 
def populate_list(file: str) -> tuple[list, int]:
    """Convert & sum phred qs / each nucleotide index(position)"""
    #intialize list(s)
    qs_list: list[float] = []
    qs_list = init_list(qs_list)
    #open FASTQ
    with gzip.open(file, mode='rt') as fh1:
        i=0
        for line in fh1:
            line = line.strip("\n")
            #isolating qs lines 
            i+=1
            if i %4==0:
                #for every position [j] in qs line
                for j,letter in enumerate(line):
                    #conv. qs letter --> qs # 
                    qscore=bioinfo.convert_phred(letter)
                    #ADD new # to sum in position j (summing q scores per position) 
                    qs_list[j]+=qscore
        return (qs_list,i)
    

my_list, num_lines = populate_list(filename)

#CALC QS AVGs
record_count=num_lines/4
#DIVISION FOR AVG 
#for every position,sum in each position in qs_list
for k,qsum in enumerate(my_list):
    qsum=qsum/record_count
    my_list[k]=qsum
    #print(f"{k}\t{qsum}")

#PLOTTING RESULTS 
import matplotlib.pyplot as plt
#setting data
avg_qs=my_list
nuc_ps=list(range(len(avg_qs)))
#setting plot 
fig, ax = plt.subplots()
#ax.plot(nuc_ps, avg_qs, 'o-', linewidth=2)
ax.bar(nuc_ps,avg_qs, log=False, color='c', width=0.5)
#adding titles 
plt.xlabel('Nucelotide Base Position')
plt.ylabel('Average Quality Score(phred_score)')
plt.title('Average Quality Scores vs. Position')
#Saving Figure 
plt.savefig(fname=f"{label}_qs_dist.png")
