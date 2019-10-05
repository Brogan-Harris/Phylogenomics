# Code to access Entrex database.
# Downloads every protein sequence (up to 100,000)

# Imports and Libraries required
import os
import sys
import time
import subprocess
from Bio import SeqIO
from Bio import Entrez
from Bio.Seq import Seq


def entrez_search_terms(E_mail, protein):

    """
    Function that returns number of protein
    terms on the entrez database
    """

    Entrez.email = E_email
    search_term = "{}".format(protein)
    handle = Entrez.esearch(db='protein', retstart = 0, retmax = 100000,
             term = search_term)
    result = Entrez.read(handle)
    count = result['Count']
    handle.close()

    return count


def entrez_search_list(search_term, no_results):

    """
    Function to search and return protein ids
    from the entrez database
    """

    handle = Entrez.esearch(db='protein', retstart = 0,
    retmax = no_results, term = search_term)
    protein_ids = Entrez.read(handle)['IdList']
    return proetin_ids

def download_genome_ids(genome_ids):

    """
    Function to download a list of genome id and write
    to a fasta file format.
    """

    for genome_id in genome_ids:
        record = Entrez.efetch(db="protein", id=genome_id, rettype="gb",
        retmode="text")
        time.sleep(0.2)

        # Write the genbank file from entrez
        filename = 'genBankRecord_{}.gb'.format(genome_id)
        print('Writing:{}'.format(filename))
        with open(filename, 'w') as f:
        	f.write(record.read()

        #Convert GenBank files to Fasta files.
        gbk_filename = 'genBankRecord_{}.gb'.format(genome_id)
        faa_filename = 'genBankRecord_{}.faa'.format(genome_id)
        input_handle  = open(gbk_filename, "r")
        output_handle = open(faa_filename, "w")

        #for seq_record in SeqIO.parse(input_handle, "genbank") :
        print "Dealing with GenBank record %s" % seq_record.id
        output_handle.write(">%s %s\n%s\n" % (seq_record.id,
                                              seq_record.description,
                                              seq_record.seq))
        output_handle.close()
        input_handle.close()
