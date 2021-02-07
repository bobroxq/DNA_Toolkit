from collections import Counter
from structures import NUCLEOTIDES, DNA_CODONS, DNA_REVERSECOMP


def validateSeq(dna_seq):
    """Checks the sequence to ensure it is a DNA string.
    Invalidates any sequence containing elements other than A, T, G, and/or C"""
    tmpseq = dna_seq.upper()
    for nuc in tmpseq:
        if nuc not in NUCLEOTIDES:
            return False
    return tmpseq


def count_nuc_frequency(dna_seq):
    """Counts the frequency of each nucleotide in a DNA sequence"""
    #  Alternatively: import collections
    #  return dict(collections.Counter(dna_seq))
    tmp_freq_dict = {"A": 0, "C": 0, "G": 0, "T": 0}
    for nuc in dna_seq:
        tmp_freq_dict[nuc] += 1
    return tmp_freq_dict


def transcription(dna_seq):
    """DNA --> RNA transcription. Replaces Thyamine with Uracil"""
    return dna_seq.replace("T", "U")


def reverse_comp(dna_seq):
    """Finds reverse compliment of DNA strand. Replaces A with T, G with C, etc"""
    #  Alternatively:
    #  mapping = str.maketrans('ATCG', 'TAGC')
    #  return dna_seq.translate(mapping)[::-1]
    return "".join([DNA_REVERSECOMP[nuc] for nuc in dna_seq])[::-1]


def gc_content(dna_seq):
    """Returns GC content of a DNA strand as a percentage"""
    return round((dna_seq.count("C") + dna_seq.count("G")) / len(dna_seq) * 100)


def gc_content_sub(dna_seq, k=20):
    """Returns GC content of a subsection of a longer strand of DNA, where k = size of each subsection"""
    res = []
    for i in range(0, len(dna_seq) - k + 1, k):
        subseq = dna_seq[i : i + k]
        res.append(gc_content(subseq))
    return res


def translate_seq(dna_seq, init_pos=0):
    """Translates nucleotide sequence into amino acid sequence"""
    return [
        DNA_CODONS[dna_seq[pos : pos + 3]]
        for pos in range(init_pos, len(dna_seq) - 2, 3)
    ]


def codon_usage(dna_seq, aminoacid):
    """Determines the frequency of a certain specified amino acid in a protein sequence"""
    tmplist = []
    for i in range(0, len(dna_seq) - 2, 3):
        if DNA_CODONS[dna_seq[i : i + 3]] == aminoacid:
            tmplist.append(dna_seq[i : i + 3])

    freq_dict = Counter(tmplist)
    total_weight = sum(freq_dict.values())
    for dna_seq in freq_dict:
        freq_dict[dna_seq] = round(freq_dict[dna_seq] / total_weight, 2)
    return freq_dict


def gen_reading_frames(dna_seq):
    """Generate the six reading frames of a DNA sequence, includng reverse compliment"""
    frames = []
    frames.append(translate_seq(dna_seq, 0))
    frames.append(translate_seq(dna_seq, 1))
    frames.append(translate_seq(dna_seq, 2))
    frames.append(translate_seq(reverse_comp(dna_seq), 0))
    frames.append(translate_seq(reverse_comp(dna_seq), 1))
    frames.append(translate_seq(reverse_comp(dna_seq), 2))
    return frames


def proteins_from_rf(aa_seq):
    """Search sequence for start and stop codons and compute all possible
    proteins in an amino acid sequence to return a list of possible proteins"""
    current_prot = []
    proteins = []
    for aa in aa_seq:
        if aa == "_":
            #  STOP accumulating amino acids if STOP codon is found
            if current_prot:
                for p in current_prot:
                    proteins.append(p)
                current_prot = []
        else:
            #  START accumulating amino acids if START codon found
            if aa == "M":
                current_prot.append("")
            for i, _ in enumerate(current_prot):
                current_prot[i] += aa

            # for i in range(len(current_prot)):
            #   current_prot[i] += aa
    return proteins


def all_proteins(dna_seq, start_pos=0, end_pos=0, ordered=False):
    """Compute all possible proteins for all possible reading frames,
    Protein search DB API can be used to pull protein info"""
    if end_pos > start_pos:
        rfs = gen_reading_frames(dna_seq[start_pos:end_pos])
    else:
        rfs = gen_reading_frames(dna_seq)

    res = []
    for rf in rfs:
        prots = proteins_from_rf(rf)
        for p in prots:
            res.append(p)

    if ordered:
        return sorted(res, key=len, reverse=True)
    return res
