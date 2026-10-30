#!/usr/bin/env python3

import getopt
import sys

def usage():
        print("Script to measure how many new bases a round added, from the alignment of the previous")
        print("and the current consensus (previous first, as for CoordinateExtraction.py)\n")
        print("The anchor is the first run of aligned previous-consensus bases on the extended side.")
        print("Bases of the current consensus before it are new, except as many as the previous consensus")
        print("also has before it: those were already in the model, only not aligned as a contiguous run.\n")
        print("Usage: ExtensionSize.py -i <Alignment> -s <Side> [Options]")
        print("\nArguments:")
        print("MANDATORY:")
        print("--input       | -i\t Alignment of previous and current consensus in fasta format")
        print("--side        | -s\t 'Left' or 'Right'")
        print("OPTIONAL:")
        print("--minRun      | -r\t Length of the run of aligned previous bases used as anchor (Default: 20)")
        print("--help        | -h\t Show this help message and exit")
        print("\nPrints the number of new bases (0 if there is no anchor or no new bases)")
        exit()

def main():
    try:
        options, remainder = getopt.getopt(sys.argv[1:], 'i:s:r:h', ['input=','side=','minRun=','help'])
    except getopt.GetoptError as err:
        print(str(err))
        usage()
        sys.exit(2)

    minRun = 20
    for opt, arg in options:
        if opt in ('-h', '--help'):
            usage()
        elif opt in ('-i', '--input'):
            Input = arg
        elif opt in ('-s', '--side'):
            Side = arg
        elif opt in ('-r', '--minRun'):
            minRun = int(arg)

    Seqs = []
    with open(Input) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                Seqs.append("")
            elif line:
                Seqs[-1] += line
    Prev, Curr = Seqs[0], Seqs[1]
    if Side == "Right":
        Prev, Curr = Prev[::-1], Curr[::-1]

    ##First column of the anchor run
    Run = 0
    Start = None
    for c, Char in enumerate(Prev):
        Run = Run + 1 if Char != "-" else 0
        if Run >= minRun:
            Start = c - minRun + 1
            break
    if Start is None:
        print(0)
        return

    PrevBefore = sum(Char != "-" for Char in Prev[:Start])
    CurrBefore = sum(Char != "-" for Char in Curr[:Start])
    print(max(0, CurrBefore - PrevBefore))

if __name__ == "__main__":
    main()
