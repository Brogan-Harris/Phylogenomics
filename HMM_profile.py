# Script in python to undertake HMM profile generation
# Uses HMMER (linux install)
# Written by Brogan Harris 18/09/2018

# Imports and Libraries required
import re
import sys
import time
import os.path
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import SearchIO

# Print dependencies
print(""" Dependencies: Biopython, HMMER (Linux install)
          Command line args: Blast database
          Format:  All files suffix *.fa"
      """)

# Parse in the blast database file name
blast_db = sys.argv[1]

# Make BLAST_db and store all records and discriptions in dictionaries.
command_3 = "makeblastdb -in " + blast_db + " -dbtype prot -out plant.db"
output = subprocess.check_output(command_3,
                                 shell=True, stderr=subprocess.STDOUT)
print("Blast database constructed")

# Dictionaries created to convert BLAST objects back into Fasta.seqs
all_seqs = {}
all_descriptions = {}
id_to_description = {}

#  Parse in BLAST database and convert records into index
BLAST_db = SeqIO.parse(blast_db, "fasta")
for record in BLAST_db:

    # Constructing Dictionary to convert BLAST HSP to fasta
    all_seqs[record.id] = str(record.seq)
    all_descriptions[record.description] = str(record.seq)
    id_to_description[record.id] = record.description

# Make HMM profile
for file_name in glob.glob("*.fa"):

    if os.path.isfile(file_name):

        # Build HMM profile
        command_1 = "hmmbuild --amino " + file_name + ".hmm " + file_name
        process = subprocess.Popen(command_1.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()

        # Emit HMM profile and write to a fasta file
        command_2 = "hmmemit " + file_name + ".hmm"
        process = subprocess.Popen(command_2.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()
        print(output)
        f = open("consensus_" + file_name, 'w' )
        f.write(output)
        f.close

        # BLAST database with emited profile and save results to xml
        try:
            blast_file = file_name + ".xml"
            command_4 = "blastp -query consensus_" + file_name + " -out " + \
                         blast_file + " -db plant.db -evalue 1e-10 -outfmt 5"
            output = subprocess.check_output(command_4,
                                             shell=True,
                                             stderr=subprocess.STDOUT)
            blast_results = SearchIO.read(blast_file, "blast-xml")

        except:
            print("error "+str(IOError))
    else:
        continue
