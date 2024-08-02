Kenny's Notebook for Dumultiplex 3 Assignment 

[on talapas...]

Data Location: /projects/bgmp/shared/2017_sequencing

TEST Input Location: /home/kenlai/bgmp/bioinfo/Bi622/Assignments/Demultiplex/TEST-input_FASTQ

Notebook/WD Location: /home/kenlai/bgmp/bioinfo/Bi622/Assignments/Demultiplex/Assignment-the-third

Creating interactive sessions: srun --account=bgmp -p bgmp -N 1 -c 4 --mem=100G -t 6:00:00 --pty bash 
_____________________
07/30/24
_____________________
* Started writing demultiplexing script: demultiplex.py 

#NOTES from class: 
* import gzip ---> gzip.open for these files. 
* leslie said to amke sure you open files separately outside of the with gzip.open block!, otherwise you don't want to open each file every singel time in the for loop 

* copying bioinfo.py --> repo 

_____________________
07/31/24
_____________________
* Writing fxns: 
```
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
```
* Tested my fxns using TEST-input_FASTQ files, worked succesfully
    * In the script I added @ end of fxns: ```#Testing my Fxns 
index_seq_to_header(R1,R2,R3,'test.fq')```
    * Ran the following in command line:```#./demultiplex.py -R1 /home/kenlai/bgmp/bioinfo/Bi622/Assignments/Demultiplex/TEST-input_FASTQ/R1.fq.gz -R2 /home/kenlai/bgmp/bioinfo/Bi622/Assignments/Demultiplex/TEST-input_FASTQ/R2.fq.gz -R3 /home/kenlai/bgmp/bioinfo/Bi622/Assignments/Demultiplex/TEST-input_FASTQ/R3.fq.gz -R4 /home/kenlai/bgmp/bioinfo/Bi622/Assignments/Demultiplex/TEST-input_FASTQ/R4.fq.gz -bars /home/kenlai/bgmp/bioinfo/Bi622/Assignments/Demultiplex/TEST-input_FASTQ/Known_Indexes.txt``` 

output of test = test.fq:
```
@test1 AAA-AAA
GATTACGT
+
9;<BA@@D
@test2 CNC-CCC
ATTCGATT
+
@@A?@?D?
@test3 AAA-TCG
TTTATGTA
+
=@@C?>=@
@test4 TCG-TCG
CATGCACC
+
12$+**&+
```

* Removed the following from the beginning and instead just did gzip when in the fxns (but still tried to b careful abt not placing inside a loop so thte file isn't opened every iteration) bc when I set it as a file object this way, with open would not accept it and watned a file "path" instead. 
```
#Opening Input Files, Setting to Variables
with gzip.open(R1g, mode='rt') as R1h,\
    gzip.open(R2g, mode='rt') as R2h,\
    gzip.open(R3g, mode='rt') as R3h,\
    gzip.open(R4g, mode='rt') as R4h:

    R1=R1h
    R2=R2h
    R3=R3h
    R4=R4h

#Do I need this step too? 
with open(known_barcodes,mode ="r") as bars:
```
* Started writing large demultiplexing function to 'actually' sort reads: 
```
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

```
* When using testfiles, executed all steps BFORE THRESHOLD FILTERING, outputs matched predicted sorting thus far (/home/kenlai/bgmp/bioinfo/Bi622/Assignments/Demultiplex/TEST-input_FASTQ)
- for R1, @test 2 (CNC-CCC) --> UNK_R1
- for R1, @test 3 (AAA_TCG) --> Hopped_R1


QUESTIONS: --> abt quality threshold 
1. do we check index qs too for filtering reads, or only that of 'actual seq.s'?
2. QS threshold setting: comparing means, medians, or per base(discard if even 1 base falls below threshold)? 
                



_____________________
08/01/24
_____________________
QUESTIONS: --> abt quality threshold 
1. do we check index qs too for filtering reads, or only that of 'actual seq.s'?
    * CHECK ONLY INDEXES, DONT NEED TO CHECKED READS BC ALIGNMENT WILL FILTER OUT POOR QSs ANYWAYS! 
2. QS threshold setting: comparing means, medians, or per base(discard if even 1 base falls below threshold)? 
    * per base, consider hamming distance. need to find a balance of a good ham_dist and qs threshold(s) 
3. Do I need to set different thresholds for R2 vs. R3? 
    * this is something you have to think about. R3 gets sequenced AFTER R2, and as a result, tends to have worse quality scores than R2. 
    * I think i will keep ham_dist same, but set R2 = 28, R3 = 26 

* Replacing this part to use a set instead of a list: 
```
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
```
```
#ISOLATING BARCODES COLUMN of known_barcodes.txt --> a set of known barcodes 
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
```
* Had to rewrite all the with open(file,'w') to write into diciontary instead to make it go faster (won't open file every time)
    --> Saved old draft with the "with open()"s as demultiplex(old).py 
    --> updated version was still called demultiplex.py 
ADDED: 
```
#Creating a Dictionary that holds the file paths (to speed up writing)
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
all_files["R1.fq"]=open("R1.fq","w")
all_files["R4.fq"]=open("R4.fq","w")
```


#Leslie said: 
* instead of making large dictionarires, you need to read through each record of all 4 files simutaneously instead. 
* can assume that the order of records is the SAME for all 4 files. 

#Need to rewrite the qs threshold lines to....
* instead of working on qs line, we need to go back and grab the qs from the R2/R3s ( or maybe create a new dictionary saying whether it passed the thresholds or not for each id )
* need to filter for BOTH R2 and R3 and different rhesholds. 

_____________________
08/02/2024
_____________________
TO DO: 
* need to filter for BOTH R2 and R3 and different rhesholds. 
* WHEN DONE, ADD ALL FXNS TO MASTER BIOINFO.py --> recopy updated version with new fxns into this repo
    * Rewrite script to call bioinfo.py functions instead of writing them all to simplify script. 

* NEW SCRIPT NAME: ```demultiplex_seq.py``` --> renamed back to ```demultiplex.py```for github convenience. 
* OLD DRAFTS OF SCRIPT: --> stored in directory ```old```
```
demultiplex_old.py #old code, didn't use dict to write, didn't open files @ same time 
demultiplex_old2.py #old code, didn't use dict to write, didn't open files @ same time 
demultiplex_old3.py #unorganized version of new script 
```
* Leslie told me I had to read all my files simutaneously instead, so I created a new script from scratch (worked w/ Claire & Varsheni): 
    * wrote the following to load files: 
    ```
    #OPENING ALL FILES 
R1=gzip.open(R1a,"rt")
R2=gzip.open(R2a,"rt")
R3=gzip.open(R3a,"rt")
R4=gzip.open(R4a,"rt")
    ```
    ```
    for line_R1,line_R2,line_R3,line_R4 in zip(R1,R2,R3,R4): #zip() connects each line-variable to its source-file
    i+=1 #incrementing record counter --> move to next record 
    line_R1=line_R1.strip() #reads same line b/w all 4 files @ a time
    line_R2=line_R2.strip() 
    line_R3=line_R3.strip()
    line_R4=line_R4.strip()
    ```

* Instead of making barcode dictionaries, I pulled my information from the simutaneously read files. 

* instead of doing with open (file) for everything, I made a dictionary where keys=file path, and values = the command to write into these files instead. This way the code takes less long to run!  (Julia helped me)
```
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
    all_files["R1.fq"]=open("R1.fq","w")
    all_files["R4.fq"]=open("R4.fq","w")

    return all_files
```
#### !Added filtering by hamming distance of 3, where R2 qs threshold = 29, R3 qs threshold = 26!
^^^ may consider changing later......

* Transfered hamming distance qs fxn & reverse comp fxn to bioinfo.py module instead to save space. 

* Created a shell script to sbatch the demultiplex.py script 

### !submitted job for first time, will check back another time to see if it worked....!
```
$ sbatch demultiplex.sh
Submitted batch job 7858367
```
--> FAILED 
* error 
    ```
    Traceback (most recent call last):
  File "/gpfs/projects/bgmp/kenlai/bioinfo/Bi622/Assignments/Demultiplex/Assignment-the-third/./demultiplex.py", line 106, in <module>
    R3_barcode = bioinfo.reverse_complement(R3_barcode_a)
                 ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/gpfs/projects/bgmp/kenlai/bioinfo/Bi622/Assignments/Demultiplex/Assignment-the-third/bioinfo.py", line 102, in reverse_complement
    rev_seq+=DNA_dict[base]
             ~~~~~~~~^^^^^^
KeyError: 'N'
    ```
* Summary: 
```
Percent of CPU this job got: 6%
Elapsed (wall clock) time (h:mm:ss or m:ss): 0:01.00
Exit status: 1
```
* SOLUTION ATTEMPT:Fixed rev_comp fxn in bioinfo.py to account for Ns: 
```
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
        'N':'N',
        'a':'t',
        't':'a',
        'c':'g',
        'g':'c',
        'n':'n',
    }
    #create empty string 
    rev_seq=""
    #complementing w/ dict. / populating string 
        #rmbr to reversed() bc the complement will want to be presented 5'->3' 
    for base in reversed(sequence):
        #append complementary base (from dict) to empty string for every base in seq. 
        rev_seq+=DNA_dict[base]
    return rev_seq

```
* Resubmitted sbatch w/ updated bioinfo.py 
```
$ sbatch demultiplex.sh
Submitted batch job 7858925
```


TO DO: 
* double check your qs threshold methods (different for R2 vs. R3, hamming distance of 3)

