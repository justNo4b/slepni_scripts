import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq





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


##  1. Открыть файл. Имя фала = Имя сиквенса
##  2. Экспорт
##
##
##
##


for gbfile in SeqIO.parse(file, "gb"):
    ## 1. Export raw seq into fasta
    basename = file.split(".")[0]
    sr = SeqRecord(Seq(gbfile.seq), id=basename, name=basename, description="")
    with open(basename + "_gb_.fasta", "w") as output_handle:
        SeqIO.write(sr, output_handle, "fasta")
    ## 2. Feature table
    with open(basename + ".ft_.txt", "w") as ft_handle:
        ft_handle.write(">Feature " + basename + "\n")
        for feature in gbfile.features:
            if (feature.type == "CDS"):
                orfStart = str(feature.location.start + 1)
                orfEnd   = str(feature.location.end)
                orfProd  = feature.qualifiers["label"][0]
                ft_handle.write(orfStart + "\t" + orfEnd + "\t" + "CDS\n")
                ft_handle.write("\t\t\t\tproduct\t" + orfProd + "\n")
                ft_handle.write("\t\t\t\ttransl_table 1\n\n")
    ## 3.Prepare source modifiers canvas
    ## Seq_ID->Collection_date->Country->Host->Strain
    with open(basename + ".sm_.txt", "w") as sm_handle:
        sm_handle.write(basename + "\t2021\tRussia: _insertREGION_\t_INSERT_HOST_\n\n")

