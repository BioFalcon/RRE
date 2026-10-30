#!/usr/bin/env python3

import getopt
import sys
import numpy as np

def usage():
        print("Script to add one round of extension to the merged model by stitching, not re-aligning.\n")
        print("Both alignments hold match columns only (one column per consensus position, as written by")
        print("HandRF.py), so the round's columns are placed onto the merged model at a known offset.\n")
        print("Usage: StitchMerge.py -m <MergedAln> -r <RoundAln> -e <Extension> -s <Side> -o <OutputPrefix>")
        print("\nArguments:")
        print("MANDATORY:")
        print("--merged      | -m\t Merged alignment so far (match columns only, fasta)")
        print("--round       | -r\t This round's alignment (match columns only, fasta)")
        print("--extension   | -e\t Number of new consensus positions at the extended edge of the round model")
        print("--side        | -s\t 'Left' (new positions are the first ones) or 'Right' (the last ones)")
        print("--output      | -o\t Output prefix (writes .aln.fa and .Consensus.fa)")
        print("OPTIONAL:")
        print("--SeqID       | -I\t Name of the consensus sequence (Default: Consensus)")
        print("--help        | -h\t Show this help message and exit")
        exit()

def ReadFasta(input):
    Names = []
    Seqs = []
    with open(input) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                Names.append(line[1:].split()[0])
                Seqs.append([])
            else:
                Seqs[-1].append(line.upper().replace(".", "-"))
    return Names, ["".join(x) for x in Seqs]

def main():
    try:
        options, remainder = getopt.getopt(sys.argv[1:], 'm:r:e:s:o:I:h', ['merged=','round=','extension=','side=','output=','SeqID=','help'])
    except getopt.GetoptError as err:
        print(str(err))
        usage()
        sys.exit(2)

    SeqID = "Consensus"
    for opt, arg in options:
        if opt in ('-h', '--help'):
            usage()
        elif opt in ('-m', '--merged'):
            Merged = arg
        elif opt in ('-r', '--round'):
            Round = arg
        elif opt in ('-e', '--extension'):
            Extension = int(arg)
        elif opt in ('-s', '--side'):
            Side = arg
        elif opt in ('-o', '--output'):
            Output = arg
        elif opt in ('-I', '--SeqID'):
            SeqID = arg

    MNames, MSeqs = ReadFasta(Merged)
    RNames, RSeqs = ReadFasta(Round)
    LenMerged = len(MSeqs[0])
    LenRound = len(RSeqs[0])
    LenNew = LenMerged + Extension

    #Offset of the round model in the new merged coordinates
    ##Left:  round positions 1..Extension are new, round position Extension+1 is merged position 1
    ##Right: round position LenRound-Extension is the last merged position
    if Side == "Left":
        MergedOffset = Extension
        RoundOffset = 0
    else:
        MergedOffset = 0
        RoundOffset = LenMerged - (LenRound - Extension)

    with open(Output + ".aln.fa", "w") as fh:
        for n, s in zip(MNames, MSeqs):
            Row = "-" * MergedOffset + s
            fh.write(">" + n + "\n" + Row + "-" * (LenNew - len(Row)) + "\n")
        for n, s in zip(RNames, RSeqs):
            #Round positions that fall before the start of the merged model are dropped
            if RoundOffset < 0:
                s = s[-RoundOffset:]
                Row = s
            else:
                Row = "-" * RoundOffset + s
            Row = Row[:LenNew]
            fh.write(">" + n + "\n" + Row + "-" * (LenNew - len(Row)) + "\n")

    #Majority-rule consensus over residues; N where there are no bases or a tie
    Rows = [("-" * MergedOffset + s).ljust(LenNew, "-") for s in MSeqs]
    for s in RSeqs:
        Row = s[-RoundOffset:] if RoundOffset < 0 else "-" * RoundOffset + s
        Rows.append(Row[:LenNew].ljust(LenNew, "-"))
    Aln = np.array([list(r) for r in Rows])
    Counts = np.stack([(Aln == b).sum(axis=0) for b in "ACGT"])
    Best = Counts.max(axis=0)
    Cons = "".join("N" if Best[j] == 0 or (Counts[:, j] == Best[j]).sum() > 1 else "ACGT"[Counts[:, j].argmax()]
                   for j in range(LenNew))
    with open(Output + ".Consensus.fa", "w") as fh:
        fh.write(">" + SeqID + "\n")
        for i in range(0, len(Cons), 80):
            fh.write(Cons[i:i + 80] + "\n")

    print("MERGED_LENGTH\t" + str(LenNew))

if __name__ == "__main__":
    main()
