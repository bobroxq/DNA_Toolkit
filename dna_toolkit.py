Nucleotides = ["A", "C", "G", "T"]

# check the sequence to make sure that it is a DNA string
def validateSeq(dna_seq):
    tmpseq = dna_seq.upper()
    for nuc in tmpseq:
        if nuc not in Nucleotides:
            return False
    return tmpseq


def countNucFrequency(dna_seq):
    tmpFreqDict = {"A": 0, "C": 0, "G": 0, "T": 0}
    for nuc in dna_seq:
        tmpFreqDict[nuc] += 1
    return tmpFreqDict
#  Alternatively: import collections
#  return dict(collections.Counter(dna_seq))