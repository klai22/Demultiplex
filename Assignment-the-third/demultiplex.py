#!/usr/bin/env python


# test_files
#./demultiplex.py -R2qs 29 -R3qs 26 -hd 3 -R1 /home/kenlai/bgmp/bioinfo/Bi622/Assignments/Demultiplex/TEST-input_FASTQ/R1.fq.gz -R2 /home/kenlai/bgmp/bioinfo/Bi622/Assignments/Demultiplex/TEST-input_FASTQ/R2.fq.gz -R3 /home/kenlai/bgmp/bioinfo/Bi622/Assignments/Demultiplex/TEST-input_FASTQ/R3.fq.gz -R4 /home/kenlai/bgmp/bioinfo/Bi622/Assignments/Demultiplex/TEST-input_FASTQ/R4.fq.gz -bars /home/kenlai/bgmp/bioinfo/Bi622/Assignments/Demultiplex/TEST-input_FASTQ/indexes.txt


# MANAGING INPUTS 
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
    parser.add_argument("-R2qs")
    parser.add_argument("-R3qs")
    parser.add_argument("-hd")
    return parser.parse_args()

#setting input args--> variables 
args=get_args()
R1a=args.R1
R2a=args.R2
R3a=args.R3
R4a=args.R4
known_barcodes=args.bars
qs_threshold_R2=int(args.R2qs)
qs_threshold_R3=int(args.R3qs)
ham_dist_threshold=int(args.hd)

# IMPORTING MODULES 
import bioinfo
import gzip

# DEFINING FXNS - script specific 
def open_files():
    '''
    Creating a dictionary that holds the file paths (instead of with open('w') bc that would take too long )
    '''
    all_files={}
    #For barcode pairs 
    for barcode in known_indexes:
        #Files for barcode-pairs (setting keys )
        name1=f"outputs/{barcode}-{barcode}_R1.fq"
        name2=f"outputs/{barcode}-{barcode}_R2.fq"
        #Opening files for writing (Setting values )
        all_files[name1]=open(name1,"w")
        all_files[name2]=open(name2,"w")
    #For UNK + HOPPED + updated R1/R4 files w/ barcoded-headers 
    all_files["outputs/UNK_R1.fq"]=open("outputs/UNK_R1.fq","w")
    all_files["outputs/UNK_R2.fq"]=open("outputs/UNK_R2.fq","w")
    all_files["outputs/Hopped_R1.fq"]=open("outputs/Hopped_R1.fq","w")
    all_files["outputs/Hopped_R2.fq"]=open("outputs/Hopped_R2.fq","w")
    return all_files
    
#OPENING ALL FILES 
R1=gzip.open(R1a,"rt")
R2=gzip.open(R2a,"rt")
R3=gzip.open(R3a,"rt")
R4=gzip.open(R4a,"rt")

#Convert known_indexes --> a set of barcodes (isolating barcodes col. of known_barcodes.txt)
#initialize empty set
known_indexes=set()
#open known indexes files 
with open(known_barcodes,'r') as file:
    #skipping header line 
    next(file)
    #looping through each line & splitting by \t 
    for line in file: 
        line_parts = line.strip().split('\t')
        #append to list (column 5 (4 bc 0-based) has barcode seq.s)
        known_indexes.add(line_parts[4])

# INITIALIZING MASTER OBJ.s 
#Creating files dict. (so we can write / categorize the reads into new files)
all_files=open_files()


# DEMULTIPLEXING 
#Initalizing record-variables 
record_R1="" #empty string to hold R1's record (4 lines @ a time) during each loop 
record_R4="" #empty string to hold R4's record during each loop 
#initalizing counter-variables 
counter = 0 #tracks which line we're on (within current record)
i =0 # tracks which record we're on ]

#READING FILES (simutaneously) --> FOR EVERY LINE in FILE.......
for line_R1,line_R2,line_R3,line_R4 in zip(R1,R2,R3,R4): #zip() connects each line-variable to its source-file
    i+=1 #incrementing record counter --> move to next record 
    line_R1=line_R1.strip() #reads same line b/w all 4 files @ a time
    line_R2=line_R2.strip()
    line_R3=line_R3.strip()
    line_R4=line_R4.strip()

    #ISOLATING INDEX BARCODES & QS LINES for current record
    #Isolate Seq. line 
    if i%4==2:
        R2_barcode = str(line_R2)
        R3_barcode_a = str(line_R3)
        #Rev. Comp. R3_barcode 
        R3_barcode = bioinfo.reverse_complement(R3_barcode_a)
    #Isolate QS line 
    if i%4==0:
        R2_qs=str(line_R2)
        R3_qs=str(line_R3)

    #ISOLATING RECORD (1 @ a time): reading lines (loop) until we've built a record (4 lines)
    #Appending current line to record_R(x)-variables until end of current record is reached (counter =4)
    if counter < 4: 
        record_R1+=f"{line_R1}\n"
        record_R4+=f"{line_R4}\n"
        counter+=1 #increment line counter --> move to next line 
        #once we each the end of current record (4th line, counter = 4), reset line count bfore incremeneting to nxt record via i+=1^^^^
    elif counter ==4:
        counter=1

        #UPDATE HEADERS: Now that we've isolated a full record, APPEND BARCODE PAIRS TO HEADERS 
        #Creating barcode-pair str 
        barcode_pair=f"{R2_barcode}-{R3_barcode}"
        #Creating updated headers (ID w/ barcode-pair) {splitting by \n, appending barcode_pair to ID line[0]}
            #Resetting Variables 
        ID1=""
        ID4=""
        R1_OD=""
        R4_ID=""
            #Updating Variables
        ID1=record_R1.split('\n')[0]
        ID4=record_R4.split('\n')[0]
        R1_ID=f"{ID1} {barcode_pair}"
        R4_ID=f"{ID4} {barcode_pair}"

        #Replacing current headers (ID only, no barcodes) OF CURRENT RECORD w/ updated headers 
            #Gathering non header lines from records to create NEW RECORD OBJ.(s) 
        R1_seq=record_R1.split('\n')[1]
        R1_plus=record_R1.split('\n')[2]
        R1_qs=record_R1.split('\n')[3]
        R4_seq=record_R4.split('\n')[1]
        R4_plus=record_R4.split('\n')[2]
        R4_qs=record_R4.split('\n')[3]
            #Create new record obj.s w/  updated headers (barcode pairs)
        record_R1_updated = f"{R1_ID}\n{R1_seq}\n{R1_plus}\n{R1_qs}"
        record_R4_updated = f"{R4_ID}\n{R4_seq}\n{R4_plus}\n{R4_qs}"

        #DEMULTIPLEX/FILTERING: sort record into appropriate file(s) 
        #Filtering: 
        #UNK_R(x) files 
            #if 3 bases in R2 < 29 or 3 bases in R3 qs < 26 OR if R2 or R3 barcodes NOT in known indexes, write record --> UNK_R(x).fq
        if bioinfo.hamdist_qs(R2_qs,ham_dist_threshold, qs_threshold_R2)==False or bioinfo.hamdist_qs(R3_qs,ham_dist_threshold, qs_threshold_R3)==False or R2_barcode not in known_indexes or R3_barcode not in known_indexes: 
            #writing into file via all_files dict.
            all_files["outputs/UNK_R1.fq"].write(f"{record_R1_updated}\n")
            all_files["outputs/UNK_R2.fq"].write(f"{record_R4_updated}\n")
        #Hopped_R(x).fq 
            #if R2/R3 barcodes DO NOT MATCH, read is considred a hopped read (mismatched indexes)
        elif str(R2_barcode) != str(R3_barcode):
            all_files["outputs/Hopped_R1.fq"].write(f"{record_R1_updated}\n")
            all_files["outputs/Hopped_R2.fq"].write(f"{record_R4_updated}\n")
            # all QC-filters passed at this point! 
        #R2-R3_R(x).fq
            #sorting records by sample(R2-R3 barcode pairs)
        else:
            all_files[f"outputs/{R2_barcode}-{R3_barcode}_R1.fq"].write(f"{record_R1_updated}\n")
            all_files[f"outputs/{R2_barcode}-{R3_barcode}_R2.fq"].write(f"{record_R4_updated}\n")
        #RESTTING RECORD & STARTING NEW RECORD
        record_R1 = ""
        record_R4 = ""
        record_R1+=f"{line_R1}\n"
        record_R4+=f"{line_R4}\n"

#PROCESSING THE LAST record, all variables arte stores, just never got to run bc no final line to push past the elif counter ==4: gate
if counter > 0:
#UPDATE HEADERS: Now that we've isolated a full record, APPEND BARCODE PAIRS TO HEADERS 
    #Creating barcode-pair str 
    barcode_pair=f"{R2_barcode}-{R3_barcode}"
    #Creating updated headers (ID w/ barcode-pair) {splitting by \n, appending barcode_pair to ID line[0]}
        #Resetting Variables 
    ID1=""
    ID4=""
    R1_OD=""
    R4_ID=""
        #Updating Variables
    ID1=record_R1.split('\n')[0]
    ID4=record_R4.split('\n')[0]
    R1_ID=f"{ID1} {barcode_pair}"
    R4_ID=f"{ID4} {barcode_pair}"

    #Replacing current headers (ID only, no barcodes) OF CURRENT RECORD w/ updated headers 
        #Gathering non header lines from records to create NEW RECORD OBJ.(s) 
    R1_seq=record_R1.split('\n')[1]
    R1_plus=record_R1.split('\n')[2]
    R1_qs=record_R1.split('\n')[3]
    R4_seq=record_R4.split('\n')[1]
    R4_plus=record_R4.split('\n')[2]
    R4_qs=record_R4.split('\n')[3]
        #Create new record obj.s w/  updated headers (barcode pairs)
    record_R1_updated = f"{R1_ID}\n{R1_seq}\n{R1_plus}\n{R1_qs}"
    record_R4_updated = f"{R4_ID}\n{R4_seq}\n{R4_plus}\n{R4_qs}"

    #DEMULTIPLEX/FILTERING: sort record into appropriate file(s) 
    #Filtering: 
    #UNK_R(x) files 
        #if 3 bases in R2 < 29 or 3 bases in R3 qs < 26 OR if R2 or R3 barcodes NOT in known indexes, write record --> UNK_R(x).fq
    if bioinfo.hamdist_qs(R2_qs,ham_dist_threshold, qs_threshold_R2)==False or bioinfo.hamdist_qs(R3_qs,ham_dist_threshold, qs_threshold_R3)==False or R2_barcode not in known_indexes or R3_barcode not in known_indexes: 
        #writing into file via all_files dict.
        all_files["outputs/UNK_R1.fq"].write(f"{record_R1_updated}\n")
        all_files["outputs/UNK_R2.fq"].write(f"{record_R4_updated}\n")
    #Hopped_R(x).fq 
        #if R2/R3 barcodes DO NOT MATCH, read is considred a hopped read (mismatched indexes)
    elif str(R2_barcode) != str(R3_barcode):
        all_files["outputs/Hopped_R1.fq"].write(f"{record_R1_updated}\n")
        all_files["outputs/Hopped_R2.fq"].write(f"{record_R4_updated}\n")
        # all QC-filters passed at this point! 
    #R2-R3_R(x).fq
        #sorting records by sample(R2-R3 barcode pairs)
    else:
        all_files[f"outputs/{R2_barcode}-{R3_barcode}_R1.fq"].write(f"{record_R1_updated}\n")
        all_files[f"outputs/{R2_barcode}-{R3_barcode}_R2.fq"].write(f"{record_R4_updated}\n")
    

#CLOSING ALL FILES: By using all_files dict. to write, files are opened & never closed. Need to close all files oepend by dict.
for files in all_files.values():
    files.close()
