# Script to convert fasta file to phylip

# Libraries and imports required
from Bio import AlignIO
from Bio import SeqIO
import argparse
import sys

# Print dependencies
print(""" Dependencies: Biopython
          Command line args: fasta file [1], phylip outfile [2]
          Format:  All files suffix *.aln, all species "Arabidopsis_thaliana"
      """)

# Command line arguments (input fasta, output phylip)
input_fasta = sys.argv[1]
output_phy = sys.argv[2]

# Functions

def record_formatter(records):

    """ Functon to make records consistent with phylip format.
        Takes records and returns a dictionary consiting of
        phylip compatible descriptions and sequences """

    # Store description and sequence
    desc_seq = {}

    # Loop through records and format for phylip.
    for record in records:
        record_description = record.description
        record_description = record_description.split("_")
        short_rec = record_description[0]
        short_rec = short_rec[0] + "_"
        short_rec = record_description[1]
        description = short_rec1 + short_rec[0:8]
        description = description.replace(".","")
        desc_seq[description] = record.seq
    return desc_seq

#  Reformat the fast file ready for conversion
with open(input_fasta, "rU") as input_handle:
    records = SeqIO.parse(input_handle, "fasta")
    desc_seq = record_formatter(records)

    ofile = open("temp.fasta", "w")
    for key, val in desc_seq.items():
        ofile.write(">" + key + "\n" + str(val) + "\n")
    ofile.close()

# Open the temporry remformatted fast file
input_handle = open("temp.fasta", "rU")
output_handle = open(output_phy, "w")

# Convert fasta file to phylip and write to directory.
alignments = AlignIO.parse(input_handle, "fasta")
AlignIO.write(alignments, output_handle, "phylip")
output_handle.close()
input_handle.close()
