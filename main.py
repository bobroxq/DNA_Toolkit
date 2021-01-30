from dna_toolkit import *
import random

#  Creating a random DNA sequence for testing
rndDNAStr = ''.join([random.choice(Nucleotides)
                    for nuc in range(50)])

DNAStr = validateSeq(rndDNAStr)
print(countNucFrequency(DNAStr))