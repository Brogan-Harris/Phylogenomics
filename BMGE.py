# Scipt to trim multiple sequence alignments
# Wrtten by Brogan Harris 01/05/2018

# Libraries to be imported
import sys
import subprocess

# Print dependencies
print(""" Dependencies: BMGE
          Command line args: Gene to be trimmed e.g. "GORK" [1]
                             BLOSUM matrix used for trimming - 30 or 70 [2]
          Format:  Alignment to be formated like "GORK_alignment.aln in dir."
      """)


# Functions
def bmge(aln_file, output_file):

    """ Trimms alignments with BMGE - requires Java """

    my_command = "java -jar ~/Desktop/Software/BMGE-1.12/BMGE.jar -i " + \
                  aln_file +  " -t AA -m BLOSUM" + BM + " -of " + output_file
    output = subprocess.check_output(my_command,
                                     shell=True, stderr=subprocess.STDOUT)
    return output

gene = sys.argv[1]
BM = sys.argv[2]
aln_file = gene + "_alignment.aln"
output_file = gene + "_trimmed_30.fasta"
output = bmge(aln_file, output_file)
print(output)
