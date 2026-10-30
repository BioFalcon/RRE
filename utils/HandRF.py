#!/usr/bin/env python3

import getopt
import sys
import numpy as np

def usage():
        print("Script to define match columns (#=GC RF) of an alignment for hmmbuild --hand,")
        print("and to trim the unsupported edge of an extension\n")
        print("Usage: HandRF.py -i <Alignment> -o <OutputPrefix> [Options]")
        print("\nArguments:")
        print("MANDATORY:")
        print("--input       | -i\t Alignment in fasta format")
        print("--output      | -o\t Output prefix (writes .stk, .aln.fa and .support.tsv)")
        print("OPTIONAL:")
        print("--minOcc      | -m\t Min. fraction of the sequences spanning a column that have a residue in it, for it to be a match column (Default: 0.5)")
        print("--minRes      | -r\t Min. number of residues in a column for it to be a match column (Default: 3)")
        print("--trimSide    | -s\t Side to trim unsupported columns from: 'Left', 'Right' or 'None' (Default: 'None')")
        print("--minCons     | -c\t Min. mean conservation (frequency of the most common base) over the window at the new edge (Default: 0.45)")
        print("--window      | -w\t Window size, in match columns, for the conservation check (Default: 10)")
        print("--minSeqs     | -n\t Min. number of residues for a match column to count as supported (Default: 5)")
        print("--trimCoverage| -t\t Trim match columns at both ends with fewer residues than this fraction of the median")
        print("                  \t residues per match column, e.g. 0.1 for a genome-wide model (Default: 0, no trimming)")
        print("--name        | -N\t Name for the #=GF ID line (Default: output prefix)")
        print("--help        | -h\t Show this help message and exit")
        exit()

def logmsg(input, output, minOcc, minRes, trimSide, minCons, window, minSeqs):
        print("_________________________________________")
        print("Running with the following arguments:")
        print("Input       : ", input)
        print("Output      : ", output)
        print("minOcc      : ", minOcc)
        print("minRes      : ", minRes)
        print("trimSide    : ", trimSide)
        print("minCons     : ", minCons)
        print("window      : ", window)
        print("minSeqs     : ", minSeqs)
        print("_________________________________________")
        print("")

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
    Seqs = ["".join(x) for x in Seqs]
    if len(set(len(x) for x in Seqs)) != 1:
        sys.exit("ERROR: sequences in " + input + " are not aligned (different lengths)")
    return Names, np.array([list(x) for x in Seqs])

def ColumnStats(Aln):
    #Residues (N counts as neither residue nor gap)
    IsBase = np.isin(Aln, list("ACGT"))
    IsRes  = Aln != "-"
    #Span of each sequence: first to last non-gap column (fragment aware)
    Cols = np.arange(Aln.shape[1])
    First = np.where(IsRes.any(axis=1), IsRes.argmax(axis=1), Aln.shape[1])
    Last  = np.where(IsRes.any(axis=1), Aln.shape[1] - 1 - IsRes[:, ::-1].argmax(axis=1), -1)
    Spanning = ((Cols[None, :] >= First[:, None]) & (Cols[None, :] <= Last[:, None])).sum(axis=0)
    NBase = IsBase.sum(axis=0)
    #Frequency of the most common base among bases
    Counts = np.stack([(Aln == b).sum(axis=0) for b in "ACGT"])
    MaxFreq = np.where(NBase > 0, Counts.max(axis=0) / np.maximum(NBase, 1), 0)
    return Spanning, NBase, MaxFreq

def main():
    try:
        options, remainder = getopt.getopt(sys.argv[1:], 'i:o:m:r:s:c:w:n:t:N:h', ['input=','output=','minOcc=','minRes=','trimSide=','minCons=','window=','minSeqs=','trimCoverage=','name=','help'])
    except getopt.GetoptError as err:
        print(str(err))
        usage()
        sys.exit(2)

    #Set default values
    minOcc   = 0.5
    minRes   = 3
    trimSide = 'None'
    minCons  = 0.45
    window   = 10
    minSeqs  = 5
    trimCov  = 0.0
    Name     = None

    #Parse options
    for opt, arg in options:
        if opt in ('-h', '--help'):
            usage()
        elif opt in ('-i', '--input'):
            Input = arg
        elif opt in ('-o', '--output'):
            Output = arg
        elif opt in ('-m', '--minOcc'):
            minOcc = float(arg)
        elif opt in ('-r', '--minRes'):
            minRes = int(arg)
        elif opt in ('-s', '--trimSide'):
            trimSide = arg
        elif opt in ('-c', '--minCons'):
            minCons = float(arg)
        elif opt in ('-w', '--window'):
            window = int(arg)
        elif opt in ('-n', '--minSeqs'):
            minSeqs = int(arg)
        elif opt in ('-t', '--trimCoverage'):
            trimCov = float(arg)
        elif opt in ('-N', '--name'):
            Name = arg

    logmsg(Input, Output, minOcc, minRes, trimSide, minCons, window, minSeqs)
    if Name is None:
        Name = Output.split("/")[-1]

    #Main code
    ##Read alignment
    Names, Aln = ReadFasta(Input)

    ##Define match columns
    Spanning, NBase, MaxFreq = ColumnStats(Aln)
    Occ = NBase / np.maximum(Spanning, 1)
    RF = (Occ >= minOcc) & (NBase >= minRes)

    ##Conservation of match columns; columns with too few residues count as unsupported
    Cons = np.where(NBase >= minSeqs, MaxFreq, 0)
    RFCols = np.where(RF)[0]
    Trimmed = 0
    if trimSide in ('Left', 'Right') and len(RFCols) > window:
        Order = RFCols if trimSide == 'Left' else RFCols[::-1]
        #Walk inwards from the extension edge until a window of match columns is conserved
        Keep = None
        for k in range(0, len(Order) - window + 1):
            if Cons[Order[k:k + window]].mean() >= minCons:
                Keep = k
                break
        if Keep is None:
            Keep = len(Order)
        #A window straddling the boundary can pass while its outer columns are still flank:
        #shave the edge back to the first individually conserved column
        while Keep < len(Order) and Cons[Order[Keep]] < minCons:
            Keep += 1
        Trimmed = Keep
        if Keep > 0:
            if trimSide == 'Left':
                ColsKeep = np.arange(Order[Keep], Aln.shape[1]) if Keep < len(Order) else np.array([], dtype=int)
            else:
                ColsKeep = np.arange(0, Order[Keep] + 1) if Keep < len(Order) else np.array([], dtype=int)
            Aln = Aln[:, ColsKeep]
            RF = RF[ColsKeep]
            Spanning, NBase, MaxFreq = Spanning[ColsKeep], NBase[ColsKeep], MaxFreq[ColsKeep]
            Cons = Cons[ColsKeep]
            #Drop sequences left without residues
            HasRes = (Aln != "-").any(axis=1)
            Names = [n for n, h in zip(Names, HasRes) if h]
            Aln = Aln[HasRes]

    ##Ends covered by only a few copies: trim from both ends while a match column has fewer residues
    ##than trimCov times the median over all match columns
    CovTrimmed = 0
    RFCols = np.where(RF)[0]
    if trimCov > 0 and len(RFCols) > 0:
        MinRes = trimCov * np.median(NBase[RFCols])
        Good = RFCols[NBase[RFCols] >= MinRes]
        if len(Good) > 0:
            ColsKeep = np.arange(Good[0], Good[-1] + 1)
            CovTrimmed = int(RF.sum() - RF[ColsKeep].sum())
            Aln = Aln[:, ColsKeep]
            RF = RF[ColsKeep]
            Spanning, NBase, MaxFreq = Spanning[ColsKeep], NBase[ColsKeep], MaxFreq[ColsKeep]
            HasRes = (Aln != "-").any(axis=1)
            Names = [n for n, h in zip(Names, HasRes) if h]
            Aln = Aln[HasRes]

    ##Write per-column support table
    with open(Output + ".support.tsv", "w") as fh:
        fh.write("column\tRF\tspanning\tbases\tmaxfreq\n")
        for j in range(Aln.shape[1]):
            fh.write("%d\t%d\t%d\t%d\t%.3f\n" % (j + 1, RF[j], Spanning[j], NBase[j], MaxFreq[j]))

    ##Write trimmed alignment in fasta format
    with open(Output + ".aln.fa", "w") as fh:
        for n, s in zip(Names, Aln):
            fh.write(">" + n + "\n" + "".join(s) + "\n")

    ##Write the match columns only (one column per consensus position, inserts dropped)
    with open(Output + ".match.aln.fa", "w") as fh:
        for n, s in zip(Names, Aln[:, RF]):
            if (s != "-").any():
                fh.write(">" + n + "\n" + "".join(s) + "\n")

    ##Write Stockholm with RF line
    Width = max([len(n) for n in Names] + [len("#=GC RF")]) + 1
    with open(Output + ".stk", "w") as fh:
        fh.write("# STOCKHOLM 1.0\n")
        fh.write("#=GF ID " + Name + "\n\n")
        for n, s in zip(Names, Aln):
            fh.write(n.ljust(Width) + "".join(s) + "\n")
        fh.write("#=GC RF".ljust(Width) + "".join("x" if r else "." for r in RF) + "\n")
        fh.write("//\n")

    print("MATCH_COLUMNS\t" + str(int(RF.sum())))
    print("TRIMMED_MATCH_COLUMNS\t" + str(Trimmed))
    print("COVERAGE_TRIMMED_MATCH_COLUMNS\t" + str(CovTrimmed))
    print("SEQUENCES\t" + str(len(Names)))

if __name__ == "__main__":
    main()
