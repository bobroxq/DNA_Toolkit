import random

from dna_toolkit import *
from utilities import colored

#  Creating a random DNA sequence for testing
rndDNAStr = ''.join([random.choice(Nucleotides)
                    for nuc in range(50)])

DNAStr = validateSeq(rndDNAStr)

print(f'\nSequence:\n{colored(DNAStr)}\n')
print(f'[1] + Sequence Length: {len(DNAStr)}\n')
print(f'[2] + Nucleotide Frequency: {countNucFrequency(DNAStr)}\n')
print(f'[3] + DNA/RNA Transcription: {transcription(DNAStr)}\n')

#  Creates diagram of double-stranded DNA
print(f"[4] + DNA String + Reverse Compliment:\n5' {DNAStr} 3'")
print(f"   {''.join(['|' for c in range(len(DNAStr))])}")
print(f"3' {reverse_comp(DNAStr)[::-1]} 5'\n")

print(reverse_comp(DNAStr))
