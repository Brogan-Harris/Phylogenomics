# Script to rename all files into an array for PBS parallel alignment.
# Written by Brogan Harris - 15/08/2019

# Imports and Libraries required #
import os
import sys

# Print dependencies
print(""" Dependencies: os
          Command line args: format to change
          Format:  All files suffix *.fasta
      """

# Store number of files in directory.
count = 0

# Read in the alignments
for alingment in glob.glob("*.fasta*"):
    new_file = str(count) + ".faa"
    os.rename(alingment, new_file)
    count += 1
