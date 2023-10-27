import argparse
from Bio import SeqIO
from Bio import SeqRecord
from Bio.Seq import Seq
from collections import namedtuple




def _fillLengthArray(filename):
    gSeqLength.clear()
    for sequence in SeqIO.parse(filename, "fasta"):
        l = len(sequence.seq)
        if (verbose): print(sequence.name + " length:" + str(l))
        gSeqLength.append(l)
    return


#################
# Calculate and present distribution info
#
def _showDistributionInfo():
    print("-------------------------------------------------------")
    ######
    # Calculate Average
    sum = 0
    i = 0
    for l in gSeqLength:
        sum = sum + l
        i = i + 1

    avg = sum / i

    print("Total sequences: " + str(i))
    print("Avg length: " + str(avg))
    print("-------------------------------------------------------")
    #################
    # Calculate median etc
    #
    gSeqLength.sort()

    num = int (len(gSeqLength) / 2)
    print("Median: " + str(gSeqLength[num]) + " Number: " + str(num) + " Compared to Avg: " + str(gSeqLength[num]/avg))
    num = int (len(gSeqLength) / 4)
    print("  0.25: " + str(gSeqLength[num]) + " Number: " + str(num) + " Compared to Avg: " + str(gSeqLength[num]/avg))
    num = int (len(gSeqLength) / 6)
    print("  0.18: " + str(gSeqLength[num]) + " Number: " + str(num) + " Compared to Avg: " + str(gSeqLength[num]/avg))
    num = int (len(gSeqLength) / 8)
    print(" 0.125: " + str(gSeqLength[num]) + " Number: " + str(num) + " Compared to Avg: " + str(gSeqLength[num]/avg))
    num = int (len(gSeqLength) / 16)
    print(" 0.062: " + str(gSeqLength[num]) + " Number: " + str(num) + " Compared to Avg: " + str(gSeqLength[num]/avg))

    print("Lowest: " + str(gSeqLength[0]) + " Number: " + str(0) + " Compared to Avg: " + str(gSeqLength[0]/avg))
    return


#################
# Step 0. Grap command line parameters
#
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--file_input", type=str)
parser.add_argument("-o", "--file_output", type=str)
parser.add_argument("-vi", "--verboseInput", action="store_true")
parser.add_argument("-vd", "--verboseDelete", action="store_true")

file        = parser.parse_args().file_input
outFile     = parser.parse_args().file_output
verbose     = parser.parse_args().verboseInput
vdelete     = parser.parse_args().verboseDelete



gSeqLength = []


#################
# Step 1. Read sequences from file, append length to the array

_fillLengthArray(file)

######################
# Step 2. Show user the data on the seq length distribution
# Currently: Average + Median and lower percentiles

_showDistributionInfo()

#######################
# Step 3. Take input about how much to cut
#
print("-------------------------------------------------------")
print("Cut sequences less than INPUT in length")
treshHold = int(input())
outFile = "higher" + str(treshHold) + "_" + outFile


####################
# Step 4. Write a new file with all short sequences cut out
#

output_handle = open(outFile, "w")


for sequence in SeqIO.parse(file, "fasta"):
    l = len(sequence.seq)
    if (l > treshHold): SeqIO.write(sequence, output_handle, "fasta")
    if (vdelete and l <= treshHold): print("Left out: " + sequence.name + " length:" + str(l))

output_handle.close()


####################
# Step 5. Analyse resulting sequence file
#
print()
print()
print("-------------------------------------------------------")
print("RESULT:")
print("-------------------------------------------------------")
_fillLengthArray(outFile)
_showDistributionInfo()


