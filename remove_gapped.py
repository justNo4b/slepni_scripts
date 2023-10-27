import argparse
from Bio import SeqIO
from Bio import AlignIO
from Bio import SeqRecord
from Bio.Seq import Seq
from collections import namedtuple


def _isSequenceUnique(sequence):
    if (len(gSeqSequence) == 0): return True

    for s in gSeqSequence:
        if (s == sequence): return False

    return True

def _calculateBadFunction(s):
    l = len(s.seq)
    gap = s.seq.count("-")
    xs  = s.seq.count("X")
    return int(((gap + xs) * 100) / l)


def _fillLengthArray(filename, mode):
    gSeqLength.clear()
    for sequence in SeqIO.parse(filename, mode):
        l = len(sequence.seq)
        gap = sequence.seq.count("-")
        xs  = sequence.seq.count("X")
        bad = _calculateBadFunction(sequence)
        if (verbose): print(sequence.name + " gaps: " + str(gap) + " X: " + str(xs) + " bad% : "+ str(bad))
        gSeqLength.append(bad)
    return



def _printWorstNumbers(startAt):
    print("-------------------------------------------------------")
    print("Worst members:")
    print("-------------------------------------------------------")
    for x in range((startAt + 2), -1, -1):
        print("bad% value: " + str(gSeqLength[x]) + " Member Num: " + str(x))


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
    print("Avg bad%: " + str(avg))
    print("-------------------------------------------------------")
    #################
    # Calculate median etc
    #
    gSeqLength.sort(reverse=True)

    num = int (len(gSeqLength) / 2)
    print("Median bad%: " + str(gSeqLength[num]) + " Member Num: " + str(num))
    num = int (len(gSeqLength) / 4)
    print("  0.25 bad%: " + str(gSeqLength[num]) + " Member Num: " + str(num))
    num = int (len(gSeqLength) / 6)
    print("  0.18 bad%: " + str(gSeqLength[num]) + " Member Num: " + str(num))
    num = int (len(gSeqLength) / 8)
    print(" 0.125 bad%: " + str(gSeqLength[num]) + " Member Num: " + str(num))
    num = int (len(gSeqLength) / 16)
    print(" 0.062 bad%: " + str(gSeqLength[num]) + " Member Num: " + str(num))

    if (num > 0): _printWorstNumbers(num)

    return


#################
# Step 0. Grap command line parameters
#
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--file_input", type=str, help="Input file.")
parser.add_argument("-o", "--file_output", type=str, help="Output file")
parser.add_argument("-vi", "--verboseInput", action="store_true", help="Print info on input sequences")
parser.add_argument("-vd", "--verboseDelete", action="store_true", help="Print info on deleted sequences")
parser.add_argument("-phylip", "--convert_to_phylip", action="store_true", help="Provide additional file in a phylip format")
parser.add_argument("-f", "--filter_same", action="store_true", help="Filter identical sequences out of alignment")

file        = parser.parse_args().file_input
outFile     = parser.parse_args().file_output
verbose     = parser.parse_args().verboseInput
vdelete     = parser.parse_args().verboseDelete
toPhylip    = parser.parse_args().convert_to_phylip
filterOn    = parser.parse_args().filter_same



gSeqLength = []
gSeqSequence = []

oMode = "fasta"

#################
# Step 1. Read sequences from file, append length to the array

_fillLengthArray(file, "fasta")

######################
# Step 2. Show user the data on the seq length distribution
# Currently: Average + Median and lower percentiles

_showDistributionInfo()

#######################
# Step 3. Take input about how much to cut
#
print("-------------------------------------------------------")
print("Cut sequences with MORE bad% than INPUT:")
treshHold = int(input())
outFile = "less_" + str(treshHold) + "_" + outFile


####################
# Step 4. Write a new file with all short sequences cut out
#

output_handle = open(outFile, "w")


for sequence in SeqIO.parse(file, "fasta"):
    bad = _calculateBadFunction(sequence)
    uniq = _isSequenceUnique(sequence.seq)
    if (vdelete and bad >= treshHold): print("Left out: " + sequence.name + " bad:" + str(bad))
    if (vdelete and uniq == False): print("Left out: " + sequence.name + " -> Non-Unique")
    if (bad < treshHold and uniq):
        SeqIO.write(sequence, output_handle, "fasta")
        gSeqSequence.append(sequence.seq)


output_handle.close()

#####
# If conversion needed, convert to phylip
#
# Usefull for phyML pipeline
if (toPhylip): AlignIO.convert(outFile, "fasta", outFile + ".phylip", "phylip-relaxed")


####################
# Step 5. Analyse resulting sequence file
#
print()
print()
print("-------------------------------------------------------")
print("RESULT:")
print("-------------------------------------------------------")
_fillLengthArray(outFile, "fasta")
_showDistributionInfo()


