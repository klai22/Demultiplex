#!/usr/bin/env python

#./demux_stats.py

#Will need to activate conda env w/ matplotlib BFORE EXECUTING! 
    #conda activate bgmp_py312

#Importing barcode pair counts table (output by demultiplex.sh)
with open("Barcode_Pair_Counts.txt","r") as fh: 
    #PERCENTAGE OF READS FOR EACH SAMPLE 
    fh.readline() #skipping header 
    #initalize a dict to hold barcode_pair [key] and percentages [value], one for samples, one for swapped indexes (hopped indexes)
    barcode_count_dict={}
    swapped_dict={}
    #initalize variable to calc. total % of swapped indexes vs. not 
    swapped_count=float(0)
    valid_samples=float(0)
    unk_count=float(0)
    #calc. sum of all read counts 
    #initalize sum obj. 
    read_count_sum=float(0)
    #looping through lines to keep count 
    for line in fh:
        line_parts=line.strip().split('\t')
        read_count=float(line_parts[2])
        read_count_sum+=read_count
    
    #populating dict.s / variables
        if line_parts[0]!=line_parts[1]:
            swapped_dict[f"{line_parts[0]}-{line_parts[1]}"]=(float(line_parts[2])/float(read_count_sum))*100
            swapped_count+=float(line_parts[2])
        elif line_parts[0]=="unknown" or line_parts[1]=="unknown":
            unk_count+=float(line_parts[2])
        elif line_parts[0]==line_parts[1]: 
            barcode_count_dict[f"{line_parts[0]}-{line_parts[1]}"]=(float(line_parts[2])/float(read_count_sum))*100
            valid_samples+=float(line_parts[2])
    
    #Final Calc.s for Overall Pie Chart 
    Perc_Sample=float((valid_samples/read_count_sum)*100)
    Perc_Hopped=float((swapped_count/read_count_sum)*100)
    Perc_Unk=float((unk_count/read_count_sum)*100)

# VISUALIZATION 
import matplotlib.pyplot as plt
#Pie Chart: % of reads from each sample (excluded hopped)
s_labels = list(barcode_count_dict.keys())
s_sizes = barcode_count_dict.values()

fig, ax = plt.subplots()
ax.pie(s_sizes, labels=s_labels,autopct='%1.1f%%')

#adding titles 
ax.set_title("Percentage of Reads from each Sample (excluding hopped)")
#Saving Figure 
plt.savefig(fname="Sample_Distribution.png")

#Pie Chart: % of reads from each hopped_index
h_labels = list(swapped_dict.keys())
h_sizes = swapped_dict.values()

fig, ax = plt.subplots()
ax.pie(h_sizes, labels=h_labels,autopct='%1.1f%%')

#adding titles 
ax.set_title("Percentage of Reads from each Hopped_Index")
#Saving Figure 
plt.savefig(fname="Swapped_Index_Distribution.png")

#Pie Chart: % of reads from overall categories
o_labels = ["Samples","Hopped_Indexes","Unknown"]
o_sizes = [Perc_Sample,Perc_Hopped,Perc_Unk] 

fig, ax = plt.subplots()
ax.pie(o_sizes, labels=o_labels,autopct='%1.1f%%')

#adding titles 
ax.set_title("Percentage Samples vs. Swapped_Indexes vs. Unknown")
#Saving Figure 
plt.savefig(fname="Overall_Distribution.png")


#Bar Plot (Horizontal): % of reads from each sample (excluded hopped)
fig, ax = plt.subplots(figsize=(14,8)) #extending width of figure bc barcodes get cutoff 
s_barcodes = list(barcode_count_dict.keys())
s_percentages = barcode_count_dict.values()
ax.barh(s_barcodes, s_percentages)
#adding titles 
ax.set_xlabel('Percentage of Reads')
ax.set_title('Percentage of Reads from each Sample (excluding hopped)')
ax.invert_yaxis()

#Saving Figure 
plt.savefig(fname="Sample_Distribution_barh.png")


#Bar Plot (Horizontal): % of reads from each hopped_index
fig, ax = plt.subplots(figsize=(14,19)) #extending width of figure bc barcodes get cutoff 
h_barcodes = list(swapped_dict.keys())
h_percentages = swapped_dict.values()
ax.barh(h_barcodes, h_percentages)
#adding titles 
ax.set_xlabel('Percentage of Reads')
ax.set_title('Percentage of Reads from each Hopped_Index')
ax.invert_yaxis()

#Saving Figure 
plt.savefig(fname="Swapped_Index_Distribution_barh.png")


#THE CODE BELOW WAS USED TO REPORT ALL SAMPLES, MIXING HOPPED AND UNHOPPED 

# #Importing barcode pair counts table (output by demultiplex.sh)
# with open("Barcode_Pair_Counts.txt","r") as fh: 
#     #PERCENTAGE OF READS FOR EACH SAMPLE 
#     fh.readline() #skipping header 
#     #initalize a dict to hold barcode_pair [key] and percentages [value] 
#     barcode_count_dict={}
#     #calc. sum of all read counts 
#     #initalize sum obj. 
#     read_count_sum=int(0)
#     #looping through lines to keep count 
#     for line in fh:
#         line_parts=line.strip().split('\t')
#         read_count=int(line_parts[2])
#         read_count_sum+=read_count
#     #populating dict 
#         barcode_count_dict[f"{line_parts[0]}-{line_parts[1]}"]=(int(line_parts[2])/int(read_count_sum))*100
    
# # VISUALIZATION 

# #Pie Chart: % of reads from each sample (including hopped)
# import matplotlib.pyplot as plt

# labels = list(barcode_count_dict.keys())
# sizes = barcode_count_dict.values()

# fig, ax = plt.subplots()
# ax.pie(sizes, labels=labels,autopct='%1.1f%%')

# #adding titles 
# ax.set_title("Percentage of Reads from each Sample (including hopped)")
# #Saving Figure 
# plt.savefig(fname="Sample Distribution.png")









