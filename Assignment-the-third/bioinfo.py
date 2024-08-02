#!/usr/bin/env python

# Author: kenneth.lok.sun.lai@gmail.com

# Check out some Python module resources:
#   - https://docs.python.org/3/tutorial/modules.html
#   - https://python101.pythonlibrary.org/chapter36_creating_modules_and_packages.html
#   - and many more: https://www.google.com/search?q=how+to+write+a+python+module

'''This module is a collection of useful bioinformatics functions
written during my KCGIP Bioinformatics and Genomics Program coursework.'''

__version__ = "0.5"         # Read way more about versioning here:
                            # https://en.wikipedia.org/wiki/Software_versioning

def convert_phred(letter: str) -> int:
    '''Converts a single character into a phred score'''
    return ord(letter) -33

def qual_score(phred_score: str) -> float:
    """calc. avg. qs for input qs line"""
    #setting init. variables
    score_sum=0
    for letter in phred_score:
        score=convert_phred(letter)
        #adding new score to initialized variable (sum_score)
        score_sum=score_sum+score
    #dividing
    return score_sum/len(phred_score)

def validate_base_seq(seq:str,RNAflag:bool=False)->bool:
    '''This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive.'''
    DNAbases = set('ATGCNatcgn')
    RNAbases = set('AUGCNaucgn')
    return set(seq)<=(RNAbases if RNAflag else DNAbases)

def gc_content(DNA):
    '''Returns GC content of a DNA or RNA sequence as a decimal between 0 and 1.'''
    assert validate_base_seq(DNA), "String contains invalid characters - are you sure you used a DNA sequence?"
    DNA = DNA.upper()
    return (DNA.count("G")+DNA.count("C"))/len(DNA)

def calc_median(lst: list) -> float:
    '''Given a sorted list, returns the median value of the list'''
    #if length of var is even %2=0
    if len(lst)%2==0:
        #a = return list[i] where i = lenth of var /2
        a = lst[(int((len(lst)/2)))-1]
        #b = return list[i] where i = ( lenth of var /2 ) + 1
        b = lst[(int(((len(lst)/2)+1)))-1]
        #return (a+ b)/2 
        return (a+b)/2
    #else (odd) 
    else:
        #return list[i] where i = lenth of var + 1 / 2
        return lst[(int(((len(lst)+1)/2)))-1]
        #Had to add -1 for each of the index #s bc when calling lst[i], it counts the 1st index as 0

#filename=fasta_file (pep.all.fa), output_file=reformatted fastafile (PS7,Bi621)
def oneline_fasta(filename,output_file):
    '''compresses all seq.lines (assocaited under same header) into one line'''
    with open(filename,"r") as fh, open(output_file,"w") as oh:
        #initializing first_time boolean, saying input lines will be "first time" until specificied as false
        first_time=True 
        for line in fh:
            #if it is the first time, add line to of, and it is so far since first_time=True so far
            if first_time: 
                oh.write(line)
            #since we added the 1st line now, we no longer want first_time=True
                first_time=False
            #if line has > (header lines), add a \n to it (so its not on same line as last seq)
            elif ">" in line:
                oh.write(f"\n{line}")
            #otherwise strip \n from seq.lines 
            else:
                stripped_line=line.strip("\n")
                oh.write(stripped_line)

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

def hamdist_qs(qs_line: str,ham_dist_threshold: int, qs_threshold: int):
    '''
    Assesses whether or not current qs_line meets qs-threshold depending on hamming_distance (# of pos.s permitted to fall below qs-threshold)
    Input: qs_line (str)[normally from R2 or R3 for demultiplex.], ham_dist threshold, qs_threshold 
    Output: a boolean (T/F) saying whether or not given qs_line passed the quality-score threshold 
    #IMPORTANT: REQUIRES bioinfo.convert_phred ( a fxn that converts a qs-line from ASCII-->actual phred scores)
    '''
    #Converting QS line (ASCII-->Phred (assumes +33 encoding))
    conv_qs=[]
    for letter in qs_line:
        phred_score=convert_phred(letter)
        conv_qs.append(phred_score)
    #Calc. Ham_Dist (a counter for # of bases w/ a qs that fell BELOW threshold)
        #Init. Ham_Dist 
    ham_dist=int(0)
        #cac. ham_dist based on qs threshold 
    for score in conv_qs:
        if score <=qs_threshold:
            ham_dist+=1
        #Filtering based on threshold
    if ham_dist >= ham_dist_threshold:
        return False
    elif ham_dist < ham_dist_threshold:
        return True

if __name__ == "__main__":
    # write tests for functions above, Leslie has already populated some tests for convert_phred
    # These tests are run when you execute this file directly (instead of importing it)
    assert convert_phred("I") == 40, "wrong phred score for 'I'"
    assert convert_phred("C") == 34, "wrong phred score for 'C'"
    assert convert_phred("2") == 17, "wrong phred score for '2'"
    assert convert_phred("@") == 31, "wrong phred score for '@'"
    assert convert_phred("$") == 3, "wrong phred score for '$'"
    assert convert_phred("E") == 36, "wrong phred score for 'E'"
    print("Your convert_phred function is working! Nice job")

    assert validate_base_seq("AATAGAT") == True, "Validate base seq does not work on DNA"
    assert validate_base_seq("AAUAGAU", True) == True, "Validate base seq does not work on RNA"
    assert validate_base_seq("Hi there!") == False, "Validate base seq fails to recognize nonDNA"
    assert validate_base_seq("Hi there!", True) == False, "Validate base seq fails to recognize nonDNA"
    assert validate_base_seq("AATAGAT"), "Validate base seq does not work on DNA"
    assert validate_base_seq("AAUAGAU", True), "Validate base seq does not work on RNA"
    assert validate_base_seq("R is the best!")==False, "Not a DNA string"
    assert validate_base_seq("aatagat"), "Validate base seq does not work on lowercase DNA"
    assert validate_base_seq("aauagau", True), "Validate base seq does not work on lowercase RNA"
    assert validate_base_seq("TTTTtttttTTT")
    assert validate_base_seq("Iliketoeatcookies")==False, "Not a DNA string"
    print("Passed DNA and RNA tests")

    assert calc_median([1,2,100]) == 2, "calc_median function does not work for odd length list"
    assert calc_median([1,2]) == 1.5, "calc_median function does not work for even length list"
    assert calc_median([1,2,3]) == 2
    assert calc_median([5,6,7,8]) == 6.5
    assert calc_median([1,1,1,1,1,1,1,1,100]) == 1
    assert calc_median([7]) == 7
    assert calc_median([50,100]) == 75
    assert calc_median([13,20,98]) == 20
    print("Median successfully calculated")

    assert gc_content("GCGCGC") == 1
    assert gc_content("AATTATA") == 0
    assert gc_content("GCATCGAT") == 0.5
    assert gc_content("CGCT") == 0.75
    print("gc_content successfully calculated")

    assert qual_score("A") == 32.0, "wrong average phred score for 'A'"
    assert qual_score("AC") == 33.0, "wrong average phred score for 'AC'"
    assert qual_score("@@##") == 16.5, "wrong average phred score for '@@##'"
    assert qual_score("EEEEAAA!") == 30.0, "wrong average phred score for 'EEEEAAA!'"
    assert qual_score("$") == 3.0, "wrong average phred score for '$'"
    assert qual_score("KEH") == 39, "wrong average phred score for 'KEH'"
    print("average quality score succesfully calculated. ")

    assert reverse_complement("GGG") == "CCC", "wrong reverse_complement for seq.'GGG'"
    assert reverse_complement("NnN") == "NnN", "wrong reverse_complement for seq.'NnN'"
    print("Your reverse_complement function is working! Nice job")

    assert hamdist_qs("@@@88",2,30) == False, "wrong ham_dist qs assesment for '@@@888'"
    print("Your hamdist_qs function is working! Nice job")