# ALE rename script
# Takes in orthogroups from orthofinder, and re-formats them for ALE
# Written by Brogan Harris 19/09/2019

# Imports and Libraries required 
import os
import sys
import glob
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import SearchIO
import pandas as pd
from collections import defaultdict


# Read in the species file
species_file = sys.argv[1]
species_dict = defaultdict(list)
with open(species_file) as f:
    for l in f:
        line = l.rstrip()
        species = line.replace(" ", "_")
        species_dict[species] = []


# Save all species.
species_set = set()
with open(species_file) as f:
        for l in f:
            species = l.rstrip()
            species_set.add(species)
print("Total number of species under analysis: ", len(species_set))

# Calulate the total number of orthogroups
total_orthogroups = 0
for orthogroup in glob.glob("*.fa"):
    total_orthogroups += 1
print(total_orthogroups)

# Loop through all the orthogroups
orthogroup_count = 0
for orthogroup in glob.glob("*.fa"):

    # Make a dictionary with all species/
    species_count = {}
    for species in species_set:
        species_count[species] = 0

    orthogroup_count += 1
    print("Orthogroup under analysis: ", orthogroup)
    print("Orthologs analysed: ", orthogroup_count, "/", total_orthogroups)
    inhandle = SeqIO.parse(orthogroup, "fasta")

    # Count the occurance of each species in the sequence file.
    descriptions_seq = {}
    species_seq = {}
    for record in inhandle:
        record_description = record.description.replace(" ", "_")
        descriptions_seq[record_description] = str(record.seq)
        for species in species_set:
            if species in record_description:
                species_seq[species] = str(record.seq)
                species_count[species] += 1


    # Reformat the sequence descriptons for ALE
    ALE_rename = {}
    for species in species_set:
        species_occurence = 0
        for description, sequence in descriptions_seq.items():
            if species in description:
                if species_count[species] > 1:
                    species_occurence += 1
                    ALE_name =  species.replace("_", ".") + "_" + \
                                str(species_occurence)
                    ALE_rename[ALE_name] = descriptions_seq[description]

    # Write the new ALE formatted sequences to a FASTA file
    file_name = "ALE_" + orthogroup
    ofile = open(file_name, "w")
    for key, val in ALE_rename.items():
        ofile.write(">" + key + "\n" + val + "\n")
    ofile.close()
