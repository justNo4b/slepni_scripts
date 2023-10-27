import argparse
from Bio import SeqIO
from Bio import SeqRecord
from Bio.Seq import Seq
from collections import namedtuple


##################
#  Function for work with sequence translation
#
def _SeqToMaxProtein(seq):
    mLenght = 0
    tmpString = seq.split('*')
    for s in tmpString:
        mLenght = max(mLenght, len(s))
    return mLenght

# define holder for our data
# it holds name, seq and its revComplement and
# orfOrient - orientation of the longest ORF
#
# if it is set to true, the current seq orientation is main ORF is on the + strand now
# otherwise it is on the - strand
#
seqData = namedtuple("seqData", "name seqOG seqRC orfOrient")
seqHolder = []
seqFinal  = []
flipCount = 0


def _isThereSame(name):
    for entry in seqHolder:
        if (entry == name):
            return True
    return False

#################
# Step 0. Grap command line parameters
#
parser = argparse.ArgumentParser()
parser.add_argument("-d", "--direction", type=str)
parser.add_argument("-i", "--file_input", type=str)
parser.add_argument("-o", "--file_output", type=str)

file        = parser.parse_args().file_input
outFile     = parser.parse_args().file_output
# see in what orientation we should align sequences
# true - main ORF in normal orientation, false - in revCompl
# use tmpOrient as simple string to bool conversion
# will result in `true` for any non-empty string
tmpOrient   = parser.parse_args().direction



#################
# Step 1. Read sequences from file, create OG and revCompl sequences
#
# For now we set orfOrient as false

for sequence in SeqIO.parse(file, "fasta"):
    seqHolder.append(sequence.name)



for sequence in SeqIO.parse(outFile, "fasta"):
    same = _isThereSame(sequence.name)
    if (same): print("Same: " + sequence.name)