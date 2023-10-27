import argparse
from Bio import SeqIO
from Bio import SeqRecord
from Bio.Seq import Seq
from collections import namedtuple



###############################
# Function to remove some tech junk

def _lineIsTechnical(l):
    if "Effective search space used" in l:
        return True
    if "Lambda      K        H" in l:
        return True
    if "Gapped" in l:
        return True
    if "Length" in l:
        return True
    return False

def _filterFasta(record, names):
    for i in range(len(names)):
        if names[i] == record.name: return True
    return False

def _containKeyword(line, words):
    for i in range(len(words)):
        if words[i] in line: return True
    return False

###################################
#  Main variables

keyWord   = ["virus", "Virus", "viria", "viridae", "phage"]
fastaName = "contigs_20.fasta"
fileName  = fastaName + "-out.txt"

results = []
names   = []
fullReport = True
totalQueries = 0
viralLikeQueries = 0


####################################
# Open file, read lines till first entry

with open(fileName, "r") as f:
    while True:
        line = f.readline()
        if not line:
            break
        if "Query=" in line:
            names.append(line.replace("Query= ", "").replace("\n", ""))
            totalQueries = totalQueries + 1
            flag = False
            q = []
            q.append(line)
            while "Effective search space used" not in line:
                line = f.readline()
                if not line: break
                ### Keyword found, mark this line and entry
                if _containKeyword(line, keyWord):
                    flag = True
                    line = line.replace("\n", "    !!!!!!!\n")
                if fullReport: q.append(line)
            ### finished, add some space and pop back non-viral name
            if not flag: del names[-1]
            if fullReport: q.append("\n\n\n")
            if flag:
                results.append(q)
                viralLikeQueries = viralLikeQueries + 1
f.close()


##################################
#   print results into file
#
with open(fileName+"_results.txt", "w") as output:
    for record in results:
        for i in range(len(record)):
            ### Skip tech info
            if not _lineIsTechnical(record[i]): output.writelines(record[i])
output.close()

#################################
#
#   filter fasta for our positive results
#

print("Total processed: "  + str(totalQueries))
print("Viral queries" + str(viralLikeQueries))

output = open(fastaName+"_filtered.fasta", "w")

for sequence in SeqIO.parse(fastaName, "fasta"):
    if _filterFasta(sequence, names): SeqIO.write(sequence, output, "fasta")

output.close()



