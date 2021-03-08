import random
import dna_toolkit.toolkit as toolkit
from dna_toolkit.utilities import colored
from dna_toolkit import structures


def main(DNAStr):

    print(f"\nSequence:\n{colored(DNAStr)}\n")
    print(f"[1] + Sequence Length: {len(DNAStr)}\n")

    print(f"[2] + Nucleotide Frequency: {toolkit.count_nuc_frequency(DNAStr)}\n")

    print(f"[3] + DNA/RNA Transcription: {toolkit.transcription(DNAStr)}\n")

    #  Creates diagram of double-stranded DNA
    print(f"[4] + DNA String + Reverse Compliment:\n5' {DNAStr} 3'")
    print(f"   {''.join(['|' for c in range(len(DNAStr))])}")
    print(f"3' {toolkit.reverse_comp(DNAStr)[::-1]} 5'\n")

    print(f"[5] + GC Content: {toolkit.gc_content(DNAStr)}%\n")

    print(
        f"[6] + GC Content in Subsection k=5: {toolkit.gc_content_sub(DNAStr, k=5)}\n"
    )

    print(f"[7] + Aminoacid Sequence from DNA: {toolkit.translate_seq(DNAStr, 0)}\n")

    print(f'[8] + Codon Frequency (L): {toolkit.codon_usage(DNAStr, "L")}\n')

    print(f"[9] + Reading Frames:")
    for frame in toolkit.gen_reading_frames(DNAStr):
        print(frame)

    print(f"\n[10] + All proteins in 6 open reading frames:")
    for prot in toolkit.all_proteins(DNAStr, 0, 0, True):
        print(f"{prot}")


if __name__ == "__main__":

    #  Creating a random DNA sequence for testing
    rndDNAStr = "".join([random.choice(structures.NUCLEOTIDES) for _ in range(50)])

    DNAStr = toolkit.validateSeq(rndDNAStr)
    main(DNAStr)
