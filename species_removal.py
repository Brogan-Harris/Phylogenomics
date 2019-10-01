# Script to remove specific species from alignments.
# Written by Brogan Harris - 15/08/2019

# Imports and Libraries required
import os
import sys
import glob
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq

# Print dependencies
print(""" Dependencies: Biopython, glob and subprocess
          Command line args: List of species
          Format:  All files suffix *.aln, all species "Arabidopsis_thaliana"
      """)

# Read in list of species, and convert to a dictionary.
species_file = sys.argv[1]
species_list = []
with open(species_file) as f:
    for l in f:
        line = l.rstrip()
        species = line.replace(" ", "_")
        species_list.append(species)

# Read in the alignments
for alingment in glob.glob("*.aln"):

    #Read in inhandle
    inhandle = SeqIO.parse(alingment, "fasta")

    # Store sequence and description
    all_seqs = {}
    all_descriptions = {}
    species_removed = {}

    # Loop through records in alignment
    for record in inhandle:

        # Add to dictionary / remove duplicates
        all_seqs[record.description] = str(record.seq)
        all_descriptions[record.description] = str(record.seq)

        # Reformat description to be the same as search key
        description = str(record.description)
        sequence = str(record.seq)
        species_removed[description] = sequence

        # Remove species from species list
        for species in species_list:
            if species in description:
                del species_removed[description]
                print("Sequence removed:", description)


    # Write new fasta files with species removed.
    file_name = "rem_" + str(alingment[:-7]) + "fasta"
    ofile = open(file_name, "w")
    for key, val in species_removed.items():
           ofile.write(">" + key + "\n" + val + "\n")
    ofile.close()

    print("\n")
