# Script to remove duplicate entries in fasta file

# Libraries and imports required
import sys
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import SearchIO

# Print dependencies
print(""" Dependencies: Biopython
          Command line args: fasta file
          Format:  Fasta file"
      """)

# Import and parse fasta file
fasta_file = sys.argv[1]
inhandle = SeqIO.parse(fasta_file, "fasta")

# Store number of sequences, and dictionary of sequences (description [Key])
count_inhandle = 0
descriptions = []
all_seqs = {}

# Loop though records and store the description[key] and sequence [value]
for record in inhandle:
    count_inhandle = count_inhandle + 1
    all_seqs[record.description] = str(record.seq)

# Store number of sequenes in the dictionary (all duplicates removed)
count_outhandle = 0
for seq in all_seqs:
    count_outhandle = count_outhandle +1

# Print number of sequences removed
print("count_inhandle: " , count_inhandle)
print("count_outhandle: " , count_outhandle)
print("Sequenced removed: " , (count_inhandle - count_outhandle))

# Write Dictionary to a new fasta file
file_name = "clear_" + fasta_file
ofile = open(file_name, "w")
for key, val in all_seqs.items():
    ofile.write(">" + key + "\n" + val +"\n")
ofile.close()
