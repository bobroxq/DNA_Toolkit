def colored(dna_seq):
    bcolors = {
        "A": "\033[92m",
        "C": "\033[94m",
        "G": "\033[93m",
        "T": "\033[91m",
        "U": "\033[91m",
        "reset": "\033[0;0m",
    }

    tmplzt = []

    for nuc in dna_seq:
        if nuc in bcolors:
            tmplzt.append(bcolors[nuc] + nuc)
        else:
            tmplzt.append(bcolors["reset"] + nuc)
    tmplzt.append(bcolors["reset"])
    return "".join(tmplzt)
