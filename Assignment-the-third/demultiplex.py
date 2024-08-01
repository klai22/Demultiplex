#!/usr/bin/env python
#BFORE RUNNING, BE SURE TO DELETE ALL EXISTING OUPUT FILES (UNK_R1.fq, hopped_R1.fq, etc.) BC IT WILL NOT OVERWRITE, IT WILL JUST APPEND TO EXISTING FILES' CONTENTS 

#./demultiplex.py -R1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz -l read1
#./part1.py -R2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz -l read2
#./part1.py -R3 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz -l index1
#./part1.py -R4 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz -l index2

#./demultiplex.py -R1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz -R2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz -R3 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz -R4 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz -bars /projects/bgmp/shared/2017_sequencing/indexes.txt

#./demultiplex.py -R1 /home/kenlai/bgmp/bioinfo/Bi622/Assignments/Demultiplex/TEST-input_FASTQ/R1.fq.gz -R2 /home/kenlai/bgmp/bioinfo/Bi622/Assignments/Demultiplex/TEST-input_FASTQ/R2.fq.gz -R3 /home/kenlai/bgmp/bioinfo/Bi622/Assignments/Demultiplex/TEST-input_FASTQ/R3.fq.gz -R4 /home/kenlai/bgmp/bioinfo/Bi622/Assignments/Demultiplex/TEST-input_FASTQ/R4.fq.gz -bars /home/kenlai/bgmp/bioinfo/Bi622/Assignments/Demultiplex/TEST-input_FASTQ/indexes.txt -qs 30

#setting get_args 
import argparse
#defining input args 
def get_args():
    parser = argparse.ArgumentParser(description="setting global vars for parsing")
    parser.add_argument("-R1")
    parser.add_argument("-R2")
    parser.add_argument("-R3")
    parser.add_argument("-R4")
    parser.add_argument("-bars")
    parser.add_argument("-qs")
    return parser.parse_args()

#setting input args--> variables 
args=get_args()
R1=args.R1
R2=args.R2
R3=args.R3
R4=args.R4
qs_threshold=int(args.qs)
known_barcodes=args.bars

# IMPORTING MODULES 
import bioinfo
import gzip

# DEFINING FXNS 
def reverse_complement(sequence: str) -> str:
    '''takes a DNA seq. & generates its complementary sequence, utilizes a dictionary including correct base pairings
    Input: DNA str (template seq.)
    Expected output: DNA str (complementary seq.)
    '''
    #creating a dict for base-pairing 
    DNA_dict ={
        'A':'T',
        'T':'A',
        'C':'G',
        'G':'C',
        'a':'t',
        't':'a',
        'c':'g',
        'g':'c',
    }
    #create empty string 
    rev_seq=""
    #complementing w/ dict. / populating string 
        #rmbr to reversed() bc the complement will want to be presented 5'->3' 
    for base in reversed(sequence):
        #append complementary base (from dict) to empty string for every base in seq. 
        rev_seq+=DNA_dict[base]

    return rev_seq

assert reverse_complement("GGG") == "CCC", "wrong reverse_complement for seq.'GGG'"
print("Your reverse_complement function is working! Nice job")

def R2_barcode_dict(barcode_fq):
    '''Creates a dict holding barcode seq.s (value) for each seq. id/ header (key)
    Input: R2 
    Output: dictionary of barcodes' records 
    '''
    with gzip.open(barcode_fq,'rt') as file:
        #initalize empty dict 
        records = {}
        #loop through records of barcodes_fq
        while True: 
            header = file.readline().strip() 
            #fxn knows to stop when out of records to process 
            if not header: 
                break 
            seq= file.readline().strip()
            #skipping + & qs lines 
            file.readline()
            file.readline()
            records[header]=seq
    return records

def R3_barcode_dict(barcode_fq):
    '''Creates a dict holding barcode seq.s (value) for each seq. id/ header (key). Accounts for need to reverse-complement R3 seq.
    Input: R3
    Output: dictionary of barcodes' records 
    '''
    with gzip.open(barcode_fq,'rt') as file:
        #initalize empty dict 
        records = {}
        #loop through records of barcodes_fq
        while True: 
            header = file.readline().strip() 
            #fxn knows to stop when out of records to process 
            if not header: 
                break 
            seq= file.readline().strip()
            rseq=reverse_complement(seq)
            #skipping + & qs lines 
            file.readline()
            file.readline()
            records[header]=rseq
    return records

def index_seq_to_header(read_fq,R2,R3,new_fq):
    '''Created new version of R1 or R4 w/ 'R2-R3' in the ID/headers of each record
    Input:R1 OR R4 [read_fq](read fq files),R2, R3 (barcode fq files)
    Expected output: name of new fq file [new_fq]
    '''
    #Create R2 & R3 record dictionaries 
    R2_records=R2_barcode_dict(R2)
    R3_records=R3_barcode_dict(R3)

    #WRITE a NEW R1/R4 file w/ updated headers(include barcodes at the end)
    #opening the read file, creating the updated version of read_file 
    with gzip.open(read_fq,'rt') as read_file, open(new_fq,'w') as output_file:
        #looping through records of read file 
        while True: 
            header=read_file.readline().strip()
            #fxn knows to stop when out of records to process 
            if not header: 
                break
            seq=read_file.readline().strip()
            plus=read_file.readline().strip()
            qs=read_file.readline().strip()
        
            #Writing the updated version of header w/ R2-R3 suffix 
                #getting barcodes associated w/ current header ID 
            R2_seq=R2_records[header]
            R3_seq=R3_records[header]
                #writing new header 
            new_header=f"{header} {R2_seq}-{R3_seq}"

            #Writing the modified record to output file 
            output_file.write(f"{new_header}\n{seq}\n{plus}\n{qs}\n")
    return output_file

#CREATING A FULL FXN TO ORGANIZE THE R1/R4 file depending on various criteria (Demultiplexing)
def demultiplex(label,read_fq,R2,R3):
    '''Sorts each read for given R1 or R4 fq (w/ barcodes in headers) into correct output files
    Input:
    R1.fq OR R4.fq [read_fq](updated read fq files from index_seq_to_header() fxn)
    R2, R3 (barcode fq files)
    label = if R1.fq --> "R1". if R4.fq --> "R2" 
    '''
    #opening the read file 
    with open(read_fq,'r') as read_file:
        #looping through records of read file 
        while True: 
            header=read_file.readline().strip()
            #fxn knows to stop when out of records to process 
            if not header: 
                break
            seq=read_file.readline().strip()
            plus=read_file.readline().strip()
            qs=read_file.readline().strip()

            #Isolating R2 & R3 barcodes from header 
                #split sequence up by SPACES
                #assumes the followin header structure: @K00337:83:HJKJNBBXX:8:1101:30086:1683 2:N:0:1 AAA-AAA
            parts = header.split()
                #isolate "R2-R3" suffix
            barcode_pair = parts[-1] #[-1] bc it is the very last part of string 
                #furthur split the barcode_pair 
            barcode_pair_parts = barcode_pair.split('-')
                #isolate R2 vs. R3 
            R2_index = barcode_pair_parts[0]
            R3_index = barcode_pair_parts[1] 

            #CATEGORIZING RECORDS/READS
            #Checking if R2 & R3 barcodes are (NOT) in known_indexes 
            if R2_index not in known_indexes and R3_index not in known_indexes:
                #appending record to UNK file ("a" appends record to file without overwritng)
                with open(f"UNK_{label}.fq","a") as unk_file: 
                    unk_file.write(f"{header}\n{seq}\n{plus}\n{qs}\n")
            else:
                #Checking if R2 barcode (doesn't)= R3 barcode 
                if R2_index != R3_index:
                    with open(f"Hopped_{label}.fq","a") as hopped_file:
                        hopped_file.write(f"{header}\n{seq}\n{plus}\n{qs}\n")
                else: 
                    #Checking if reads fall below QS threshold 
                        #Converting QS line (ASCII-->Phred (assumes +33 encoding))
                    phred_score=bioinfo.convert_phred(qs)
                        #Calc. AVG OR MEDIAN OR ___some per base method, would need to reconfig. code...
                    #final_score=bioinfo.qual_score(phred_score) OR bioinfo.calc_median(#turnphred_score-->a list) or ???
                        #Filtering based on threshold
                    if final_score < qs_threshold:
                        with open(f"UNK_{label}.fq","a") as unk_file: 
                            unk_file.write(f"{header}\n{seq}\n{plus}\n{qs}\n")
                    else: #if => qs_threshold
                        #Passed all filters, add record to unique barcode-pair labeled .fq files
                        with open(f"{barcode_pair}_{label}.fq","a") as unk_file: 
                            unk_file.write(f"{header}\n{seq}\n{plus}\n{qs}\n")



                    #QUESTIONS: 
                    #do we check index qs too for filtering reads, or only that of 'actual seq.s'?
                    #QS threshold setting: comparing means, medians, or per base(discard if even 1 base falls below threshold)? 
                
# EXECUTING FUNCTIONS 

#WRITING NEW R1 & R4 files w/ updated headers: 
index_seq_to_header(R1,R2,R3,'R1.fq')
index_seq_to_header(R4,R2,R3,'R4.fq')

#ISOLATING BARCODES COLUMN of known_barcodes.txt --> a list of known barcodes 
#initialize empty list
known_indexes=[]
#open known indexes files 
with open(known_barcodes,'r') as file:
    #skipping header line 
    next(file)
    #looping through each line & splitting by \t 
    for line in file: 
        line_parts = line.strip().split('\t')
        #append to list (column 5 (4 bc 0-based) has barcode seq.s)
        known_indexes.append(line_parts[4])


#Testing my Fxns 
index_seq_to_header(R1,R2,R3,'test.fq')
demultiplex("R1","test.fq",R2,R3)
#print(known_indexes)






#Opening Input Files 
#[with open] input .fq files: R1 & R4 (reads), R2 & R3 (indexes),known_barcodes.txt(lists known barcodes):
    #Looping through R3 barcodes & collecting associated R2 barcodes (same ID)
    #[For] every seq.line in R3:
        #R3C_variable = reverse_complement(R3 seq.)
        #R2_variable = isolated seq. from same read as current R3C _variable (taken from R2) 
    #Creating new versions of R1/R4 with IDs that include 'R2-R3 barcodes' in name
        #R1_new = index_seq_to_header(R1,R3C_variable,R2_variable)
        #R4_new = index_seq_to_header(R4,R3C_variable,R2_variable)
    #Looping through every record of updated R1/R4 & checking for matches in known_barcodes.txt
    #REPEAT ALL STEPS BELOW FOR R4 AS WELL! 
    #[For] every record in R1_new:
        #Set record equal to a variable 
            #record = current record
        #isolate R2 & R3 barcodes from ID --> check to see if they have matches in known_barcodes.txt
        #False --> write_record_to_file(record,UNK_R1.fq)
        #True -->:
            #Is R2 barcode == R3 barcode?:
                #False --> write_record_to_file(record,Hopping_R1.fq) #indicative of index hopping
                #True --> (barcodes are complementary):
                    #Isolate QS line associated w/ current header/ID
                #Calc. avg/median QS for current record 
                    #qs = calc_qs()
                    #Test if current record's QS meets threshold 
                        #Falls below threshold --> write_record_to_file(record,UNK_R1.fq)
                    #printing record in appropriate index-pair file
                        #Meets QS threshold --> write_record_to_file(record,{R2}-{R3}_R1.fq)



#NOTES: 
# import gzip ---> gzip.open for these files. 
# leslie said to amke sure you open files separately from the for loop, toehrwise you don't want to open each file every singel time in the for loop 