#!/usr/bin/env python3

import getopt
import sys
import numpy as np

def usage():
        print("Script to add one round of extension to the merged model by stitching, not re-aligning the copies.\n")
        print("The consensus of the round is aligned to the edge of the merged consensus. Positions matched in that")
        print("alignment share a column; positions present in only one of the two get a column of their own, with")
        print("gaps in the rows of the other block. Round positions before the first matched one are the extension.\n")
        print("Usage: StitchMerge.py -m <MergedAln> -r <RoundAln> -e <Extension> -s <Side> -o <OutputPrefix>")
        print("\nArguments:")
        print("MANDATORY:")
        print("--merged      | -m\t Merged alignment so far (fasta)")
        print("--round       | -r\t This round's alignment (match columns only, fasta)")
        print("--extension   | -e\t Extension measured for this round (only used to warn if the stitch disagrees)")
        print("--side        | -s\t 'Left' (new positions are the first ones) or 'Right' (the last ones)")
        print("--output      | -o\t Output prefix (writes .aln.fa, .Consensus.fa and .Columns.txt)")
        print("OPTIONAL:")
        print("--SeqID       | -I\t Name of the consensus sequence (Default: Consensus)")
        print("--help        | -h\t Show this help message and exit")
        print("\nThe call of every column ('-' for columns that are not consensus positions) is written to")
        print("<prefix>.Columns.txt and, if found next to the merged alignment, read back in the next round,")
        print("so a position called from one round is not dropped later for being gapped in the other rounds.")
        exit()

#Scores for the consensus to consensus alignment
MATCH, MISMATCH, GAPOPEN, GAPEXT = 2, -3, -6, -1

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

def Consensus(Seqs, minOcc=0.5):
    #Majority-rule consensus over residues, one character per column: N where there are no bases or a tie,
    #'-' where fewer than minOcc of the rows spanning the column (from their first to last base) have a base
    Aln = np.array([list(s) for s in Seqs])
    IsBase = np.isin(Aln, list("ACGTN"))
    Ncol = Aln.shape[1]
    First = np.where(IsBase.any(axis=1), IsBase.argmax(axis=1), Ncol)
    Last = np.where(IsBase.any(axis=1), Ncol - 1 - IsBase[:, ::-1].argmax(axis=1), -1)
    Cols = np.arange(Ncol)
    Spanning = ((First[:, None] <= Cols) & (Last[:, None] >= Cols)).sum(axis=0)
    Bases = IsBase.sum(axis=0)
    Counts = np.stack([(Aln == b).sum(axis=0) for b in "ACGT"])
    Best = Counts.max(axis=0)
    Out = []
    for j in range(Ncol):
        if Bases[j] == 0 or Bases[j] < minOcc * Spanning[j]:
            Out.append("-")
        elif Best[j] == 0 or (Counts[:, j] == Best[j]).sum() > 1:
            Out.append("N")
        else:
            Out.append("ACGT"[Counts[:, j].argmax()])
    return "".join(Out)

def Supported(Cons):
    #Consensus without the unsupported columns, and the column of each of its positions
    Idx = [j for j, c in enumerate(Cons) if c != "-"]
    return "".join(Cons[j] for j in Idx), Idx

def AlignEdge(a, b):
    #Affine-gap alignment of the round consensus a to the start of the merged consensus b.
    #Free: round positions before the alignment (the extension), merged positions before and after it.
    #Not free: round positions after it, the round must end inside the merged model.
    #Returns the alignment as a list of (round position or None, merged position or None)
    n, m = len(a), len(b)
    NEG = -1e9
    M = np.full((n + 1, m + 1), NEG)
    X = np.full((n + 1, m + 1), NEG)   #round residue against a gap
    Y = np.full((n + 1, m + 1), NEG)   #merged residue against a gap
    M[0, 0] = 0
    X[1:, 0] = 0
    Y[0, 1:] = 0
    for i in range(1, n + 1):
        ai = a[i - 1]
        for j in range(1, m + 1):
            bj = b[j - 1]
            s = 0 if ai == "N" or bj == "N" else (MATCH if ai == bj else MISMATCH)
            M[i, j] = s + max(M[i - 1, j - 1], X[i - 1, j - 1], Y[i - 1, j - 1])
            X[i, j] = max(M[i - 1, j] + GAPOPEN, X[i - 1, j] + GAPEXT, Y[i - 1, j] + GAPOPEN)
            Y[i, j] = max(M[i, j - 1] + GAPOPEN, Y[i, j - 1] + GAPEXT, X[i, j - 1] + GAPOPEN)

    #Best end in the last round position
    Ends = np.stack([M[n], X[n], Y[n]])
    State, j = np.unravel_index(Ends.argmax(), Ends.shape)
    i = n
    Pairs = [(None, q) for q in range(m - 1, j - 1, -1)]
    while i > 0 and j > 0:
        if State == 0:
            s = 0 if a[i - 1] == "N" or b[j - 1] == "N" else (MATCH if a[i - 1] == b[j - 1] else MISMATCH)
            Prev = [M[i - 1, j - 1], X[i - 1, j - 1], Y[i - 1, j - 1]]
            Pairs.append((i - 1, j - 1))
            State = int(np.argmax(Prev))
            i, j = i - 1, j - 1
        elif State == 1:
            Prev = [M[i - 1, j] + GAPOPEN, X[i - 1, j] + GAPEXT, Y[i - 1, j] + GAPOPEN]
            Pairs.append((i - 1, None))
            State = int(np.argmax(Prev))
            i -= 1
        else:
            Prev = [M[i, j - 1] + GAPOPEN, X[i, j - 1] + GAPOPEN, Y[i, j - 1] + GAPEXT]
            Pairs.append((None, j - 1))
            State = int(np.argmax(Prev))
            j -= 1
    #Round positions before the alignment are new; merged positions before it are free
    Pairs += [(p, None) for p in range(i - 1, -1, -1)]
    Pairs += [(None, q) for q in range(j - 1, -1, -1)]
    return Pairs[::-1]

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

    #Work as if extending to the left; the right side is the mirror image
    if Side == "Right":
        MSeqs = [s[::-1] for s in MSeqs]
        RSeqs = [s[::-1] for s in RSeqs]

    ColumnsFile = Merged[:-len(".aln.fa")] + ".Columns.txt" if Merged.endswith(".aln.fa") else Merged + ".Columns.txt"
    try:
        MCalls = open(ColumnsFile).read().strip()
        if len(MCalls) != len(MSeqs[0]):
            raise ValueError
        if Side == "Right":
            MCalls = MCalls[::-1]
    except (OSError, ValueError):
        MCalls = Consensus(MSeqs)
    MCons, MPos = Supported(MCalls)
    RCons, RPos = Supported(Consensus(RSeqs))

    #Only the edge of the merged model can overlap the round
    Window = min(len(MCons), len(RCons) + 300)
    Pairs = AlignEdge(RCons, MCons[:Window])
    Pairs += [(None, q) for q in range(Window, len(MCons))]
    Pairs = [(None if r is None else RPos[r], None if q is None else MPos[q]) for r, q in Pairs]

    #Report, in consensus positions
    Stitched = next(k for k, (r, q) in enumerate(Pairs) if q is not None)
    LastRound = max(k for k, (r, q) in enumerate(Pairs) if r is not None)
    NewPositions = sum(r is not None for r, q in Pairs[:Stitched])
    Matched = sum(r is not None and q is not None for r, q in Pairs)
    GapsInMerged = sum(r is not None and q is None for r, q in Pairs[Stitched:LastRound + 1])
    GapsInRound = sum(r is None and q is not None for r, q in Pairs[Stitched:LastRound + 1])

    #Unsupported columns of either block keep their place, between their neighbours
    Columns = []
    NextR, NextM = 0, 0
    for r, q in Pairs:
        if r is not None:
            Columns += [(x, None) for x in range(NextR, r)]
            NextR = r + 1
        if q is not None:
            Columns += [(None, x) for x in range(NextM, q)]
            NextM = q + 1
        Columns.append((r, q))
    Columns += [(x, None) for x in range(NextR, len(RSeqs[0]))]
    Columns += [(None, x) for x in range(NextM, len(MSeqs[0]))]

    MIdx = np.array([-1 if q is None else q for r, q in Columns])
    RIdx = np.array([-1 if r is None else r for r, q in Columns])

    def Place(s, Idx):
        Row = np.array(list(s) + ["-"])
        return "".join(Row[Idx])

    Rows = [Place(s, MIdx) for s in MSeqs] + [Place(s, RIdx) for s in RSeqs]

    #A column present in both blocks is called from all rows; one present in only one block keeps that
    #block's call: the other block has only gaps there because its alignment holds match columns only
    NM = len(MSeqs)
    CAll, CR = Consensus(Rows), Consensus(Rows[NM:])
    Calls = "".join(CAll[j] if r is not None and q is not None else (CR[j] if r is not None else MCalls[q])
                    for j, (r, q) in enumerate(Columns))

    if Side == "Right":
        Rows = [s[::-1] for s in Rows]
        Calls = Calls[::-1]

    with open(Output + ".aln.fa", "w") as fh:
        for n, s in zip(MNames + RNames, Rows):
            fh.write(">" + n + "\n" + s + "\n")

    with open(Output + ".Columns.txt", "w") as fh:
        fh.write(Calls + "\n")

    Cons = Calls.replace("-", "")
    with open(Output + ".Consensus.fa", "w") as fh:
        fh.write(">" + SeqID + "\n")
        for i in range(0, len(Cons), 80):
            fh.write(Cons[i:i + 80] + "\n")

    print("MERGED_LENGTH\t" + str(len(Cons)))
    print("STITCH_NEW_POSITIONS\t" + str(NewPositions))
    print("STITCH_MATCHED\t" + str(Matched))
    print("STITCH_GAPS_IN_MERGED\t" + str(GapsInMerged))
    print("STITCH_GAPS_IN_ROUND\t" + str(GapsInRound))
    if abs(NewPositions - Extension) > 10:
        print("WARNING: the stitch places %d new positions, the measured extension was %d" % (NewPositions, Extension))

if __name__ == "__main__":
    main()
