#!/usr/bin/env python3

import getopt
import sys
import numpy as np

def usage():
        print("Script to estimate how similar the copies of a repeat are to each other, and from it the identity")
        print("above which two copies are too similar to be independent insertions (e.g. segmental duplications)\n")
        print("Usage: PairwiseIdentity.py -i <Alignment> [Options]")
        print("\nArguments:")
        print("MANDATORY:")
        print("--input       | -i\t Alignment in fasta format")
        print("OPTIONAL:")
        print("--minShared   | -S\t Min. number of columns where both copies have a base to compare them (Default: 50)")
        print("--fraction    | -f\t Threshold = median + fraction * (1 - median) (Default: 0.5)")
        print("--min         | -m\t Lowest threshold returned (Default: 0.8, the lowest cd-hit-est accepts)")
        print("--max         | -M\t Highest threshold returned (Default: 0.98)")
        print("--help        | -h\t Show this help message and exit")
        print("\nPrints MEDIAN_IDENTITY and THRESHOLD")
        exit()

def main():
    try:
        options, remainder = getopt.getopt(sys.argv[1:], 'i:S:f:m:M:h', ['input=','minShared=','fraction=','min=','max=','help'])
    except getopt.GetoptError as err:
        print(str(err))
        usage()
        sys.exit(2)

    minShared = 50
    Fraction  = 0.5
    Low       = 0.8
    High      = 0.98
    for opt, arg in options:
        if opt in ('-h', '--help'):
            usage()
        elif opt in ('-i', '--input'):
            Input = arg
        elif opt in ('-S', '--minShared'):
            minShared = int(arg)
        elif opt in ('-f', '--fraction'):
            Fraction = float(arg)
        elif opt in ('-m', '--min'):
            Low = float(arg)
        elif opt in ('-M', '--max'):
            High = float(arg)

    Seqs = []
    with open(Input) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                Seqs.append([])
            elif line:
                Seqs[-1].append(line.upper())
    Aln = np.array([list("".join(x)) for x in Seqs])
    IsBase = np.isin(Aln, list("ACGT"))

    Ident = []
    for a in range(len(Aln) - 1):
        Shared = IsBase[a][None, :] & IsBase[a + 1:]
        NShared = Shared.sum(axis=1)
        Same = ((Aln[a][None, :] == Aln[a + 1:]) & Shared).sum(axis=1)
        Ident += list((Same / np.maximum(NShared, 1))[NShared >= minShared])

    if len(Ident) == 0:
        print("MEDIAN_IDENTITY\tNA")
        print("THRESHOLD\t" + str(High))
        return
    Median = float(np.median(Ident))
    print("MEDIAN_IDENTITY\t%.3f" % Median)
    print("THRESHOLD\t%.3f" % min(High, max(Low, Median + Fraction * (1 - Median))))

if __name__ == "__main__":
    main()
