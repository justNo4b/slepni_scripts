import argparse
import os
from Bio import SeqIO
from Bio import Entrez
from Bio import Phylo



def _filterBootstrap():
    # ask user for a Bootstraps params performed
    print("-------------------------------------------------------")
    print("Removing unreliable bootstraps:")
    print("-------------------------------------------------------")
    print("Max bootstraps performed:")
    bootstrapMax = int(input())
    print("Remove less than (%):")
    bootstrapFilterPercent = int(input())

    # iterate all nonterminal nodes
    nonTerminal = tree.find_clades(terminal=False)
    while True:
        clade = next(nonTerminal, False)
        if (not clade): break
        c = clade.confidence
        if ((c != None) and (c * 100 / bootstrapMax) < bootstrapFilterPercent):
            clade.confidence = 0
        if ((c != None) and (c * 100 / bootstrapMax) >= bootstrapFilterPercent):
            clade.confidence = bootstrapMax


    return



def _formatName(id, name, strain):
    s = "[" + id +"] " + name + " " + strain
    return s

#################
# Step 0. Grap command line parameters
#
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--file_input", type=str, help="Input file.")
parser.add_argument("-o", "--file_output", type=str, help="Output file")
parser.add_argument("-e", "--email", type=str, help="Email, required to access NCBI via Enterez")
parser.add_argument("-fb", "--filter_bootstrap", action="store_true", help="Set bootstraps lower than certain percent value to zero")
parser.add_argument("-st", "--save_temp", action="store_true", help="Use this option in order to save temp file")
parser.add_argument("-si", "--save_additional_info", action="store_true", help="Use this option to save additional info")

file        = parser.parse_args().file_input
outFile     = parser.parse_args().file_output
filterBoot  = parser.parse_args().filter_bootstrap
email       = parser.parse_args().email
saveTemp    = parser.parse_args().save_temp
saveInfo    = parser.parse_args().save_additional_info




Entrez.email = email
fetchNum = ""
tempFile = "temp.gb"
addInfo  = []


####################################
#
#   Step 1. Read input file (phylogenetic tree)
#   Get all IDs
#
#   use get_terminals() to get all terminal nodes as a list it seems (usefull for getting names)
#
#   use get_nonterminals() to get all internal nodes (NOTE: also root)


tree = Phylo.read(file, "newick")
terminals = tree.get_terminals()
for entry in terminals:
    fetchNum = fetchNum + entry.name + ","

#######################################
#
#   efetch - ids must be separated by ','
#
#  To get full record
#  for nucleotide use: rettype="gb", retmode="text"
#  for protein use:    rettype="gp", retmode="text"
# NOTE: for some reason nucleotide queries more sensitive to wrong-formatted IDs
# will need some fixes for later filtering
#
#   Step 2. Get full information from GenBank in requested gb format
handle = Entrez.efetch(db="protein", id=fetchNum, rettype="gp", retmode="text")
with open(tempFile, "w") as output:
    output.write(handle.read())





###########################################
#   Step 3: Annotate the tree names
#
#   smh record.annotations["organism"] sometimes produce additional
#   higher taxon shit, thus use features -> qualifiers
with open(tempFile) as temp:
    for record in SeqIO.parse(temp, "gb"):
        clade = next(tree.find_clades(name=record.id))
        add = ""
        if (("starin" in record.features[0].qualifiers)): add = record.features[0].qualifiers["strain"][0]
        if (("isolate" in record.features[0].qualifiers)): add = record.features[0].qualifiers["isolate"][0]
        clade.name = _formatName(record.id, record.features[0].qualifiers["organism"][0], add)
        if (saveInfo and ("host" in record.features[0].qualifiers)): addInfo.append(clade.name + "\t" + record.features[0].qualifiers["host"][0] + "\n")




###########################################
#   Step 4: Remove bootstraps if they are lower than
#           certain percentage
#
if (filterBoot): _filterBootstrap()



###########################################
#   Step 4: Save it in a proper format
#
Phylo.write(tree, outFile, "newick")

###########################################
#   Step 5: Write additional info
#
if (saveInfo):
    with open(outFile+"_info.txt", "w") as output:
        for record in addInfo:
                output.writelines(record)
    output.close()

###########################################
#   Step 6: Remove temp garbage
#
if (not saveTemp): os.remove(tempFile)