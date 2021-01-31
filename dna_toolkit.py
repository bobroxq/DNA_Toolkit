from structures import *

def validateSeq(dna_seq):
    """Checks the sequence to ensure it is a DNA string. Invalidates any sequence containing elements other than A, T, G, and/or C"""
    tmpseq = dna_seq.upper()
    for nuc in tmpseq:
        if nuc not in Nucleotides:
            return False
    return tmpseq


def countNucFrequency(dna_seq):
    """Counts the frequency of each nucleotide in a DNA sequence"""
    tmpFreqDict = {"A": 0, "C": 0, "G": 0, "T": 0}
    for nuc in dna_seq:
        tmpFreqDict[nuc] += 1
    return tmpFreqDict
#  Alternatively: import collections
#  return dict(collections.Counter(dna_seq))


def transcription(dna_seq):
    """DNA --> RNA transcription. Replaces Thyamine with Uracil"""
    return dna_seq.replace("T", "U")


def reverse_comp(dna_seq):
    """Finds reverse compliment of DNA strand. Replaces A with T, G with C, etc"""
    return ''.join([DNA_ReverseComp[nuc] for nuc in dna_seq])[::-1]