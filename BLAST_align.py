# Script to undertake BLAST search, align homologues and create Directory
# Written by Brogan Harris 01/09/2018
# Please acknowledge if used.

# Imports and libraries required
import os
import re
import sys
import time
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import SearchIO

# Print dependencies
print(""" Dependencies: Biopython, NCBI BLAST, Mafft 0.9 and BMGE
          Command line arguments: Gene list (1), BLAST database (2)
          Format: Gene + "_prot.fasta" i.e GORK_prot.fasta
      """)

# Command line arguments 1) Gene to BLAST 2) BLAST database
gene = sys.argv[1]
blast_db = sys.argv[2]

# Make BLAST database via NCBI
my_command = "makeblastdb -in " + blast_db + " -dbtype prot -out plant.db"
output = subprocess.check_output(my_command,
                                 shell=True, stderr=subprocess.STDOUT)
print "Blast database constructed"

# Loop through list of genes and BLAST against DB
for gene in gene_list:
    # Undertake BLAST search
    try:
        my_command = "blastp -query " + gene + "_prot.fasta -out " + gene +  \
                     "_BLAST.xml -db plant.db -evalue 0.001 -outfmt 5"
        output = subprocess.check_output(my_command, shell=True,
                 stderr=subprocess.STDOUT)
        print("Blast search completed for gene: " + gene)

        # Turn BLAST.db into a database which can be indexed, covert to dict.
        all_seqs = {}
        all_descriptions = {}
        id_to_description = {}
        inhandle = SeqIO.parse(blast_db, "fasta")
        for record in inhandle:
            all_seqs[record.id] = str(record.seq) #removes duplicates
            all_descriptions[record.description] = str(record.seq)
            id_to_description[record.id] = record.description

        # Open BLAST_XML
        blast_file = gene + "_BLAST.xml"
        blast_qresult = SearchIO.read(blast_file, "blast-xml")

        # Get BLAST hit ID's to retrieve sequences
        hits = []
        length_blast = len(blast_qresult)
        for i in range(length_blast):
            blast_hsp = blast_qresult[i][0]
            hits.append(blast_hsp.hit_id)

        # Write description and sequence to a fasta file to be used later.
        ofile = open(gene + "_BLAST_results.faa", "w")
        for hit in hits:
            if "BL_ORD_ID" in hit:
                continue
            else:
                ofile.write(">" +id_to_description[hit]+ "\n" + all_seqs[hit]+"\n")
        ofile.close()

    except:
        "Print gene not availble"

    # Aligning sequence using MAFFT
    answer = input("Would you like to align the BLAST hits - Y/N? ")
    if answer == "Y":
        print("Sequences are being aligned by MAFFT-settings: maxiterate 1000")
        my_command = "mafft --auto "+ gene +  \
                      "_BLAST_results.faa > " + gene + "_alignment.aln"
        output = subprocess.check_output(my_command, shell=True,
                 stderr=subprocess.STDOUT)
        print("\nSequences have been aligned:"+ gene +"_alignment.aln is ready")

        #global gene
        os.mkdir(gene + "._low_eval")
        try:
            my_command = "mv *" + gene +"_* " + gene + "._low_eval/"
            output  = subprocess.check_output(my_command, shell=True,
                      stderr=subprocess.STDOUT)
        except("ProcessError"):
            print("Files have been moved to folder under the the Queries name")


    # Make Directory to save the output files to.
    if answer == "N":
        #global gene
        os.mkdir(gene + "._low_eval")
        my_command = "mv *" + gene +"_* " + gene + "._low_eval/"
        output  = subprocess.check_output(my_command,
                  shell=True, stderr=subprocess.STDOUT)
        print("Files have been moved to folder under the the Queries name")
