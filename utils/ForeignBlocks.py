#!/usr/bin/env python3

import getopt
import random
import subprocess
import sys
import tempfile
import numpy as np

def usage():
        print("Script to find blocks of a model that belong to a younger, more abundant family picked up during the")
        print("extension (e.g. the model ran into a high copy repeat next to some of its copies).\n")
        print("A block is foreign when it is both abundant (many more genomic hits cover it than the rest of the model)")
        print("and young (its copies are much more similar to each other than the copies of the rest of the model).")
        print("Abundance alone is not enough: some regions of an ancient repeat are better conserved and covered than")
        print("others. The model is cut at the lowest-depth point between the block holding the anchor and the")
        print("nearest foreign block.\n")
        print("Usage: ForeignBlocks.py -t <nhmmer tblout> -L <model length> -g <genome.fa> [Options]")
        print("\nArguments:")
        print("MANDATORY:")
        print("--tbl         | -t\t nhmmer --tblout of the model against the genome")
        print("--length      | -L\t Model length (match states)")
        print("--genome      | -g\t Genome fasta (indexed), to extract sampled copies")
        print("OPTIONAL:")
        print("--anchor      | -a\t Model position that must be kept, e.g. the seed region (Default: centre of the model)")
        print("--maxE        | -E\t Max. E-value of the hits used (Default: 1e-5)")
        print("--window      | -w\t Window size in model positions (Default: 50)")
        print("--depthFold   | -f\t Min. depth relative to the model's median depth for a window to be abundant (Default: 3)")
        print("--identityGap | -d\t Min. identity above the model's baseline identity for a window to be young (Default: 0.15)")
        print("--sample      | -s\t Copies sampled per window to measure identity (Default: 40)")
        print("--cpu         | -c\t Threads for mafft (Default: 1)")
        print("--help        | -h\t Show this help message and exit")
        print("\nPrints the windows tested and KEEP <start> <end> (1-based model positions to keep)")
        exit()

def ReadHits(Tbl, MaxE):
    H = []
    for line in open(Tbl):
        if line.startswith("#"):
            continue
        x = line.split()
        if float(x[12]) > MaxE:
            continue
        a, b = int(x[6]), int(x[7])
        H.append((x[0], min(a, b) - 1, max(a, b), x[11], int(x[4]), int(x[5])))
    return H

def MedianIdentity(Seqs, CPU):
    #Median pairwise identity over shared bases of the copies, after aligning them with mafft
    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as fh:
        for k, s in enumerate(Seqs):
            fh.write(">s%d\n%s\n" % (k, s))
        Name = fh.name
    Out = subprocess.run(["mafft", "--quiet", "--thread", str(CPU), "--auto", Name], capture_output=True, text=True).stdout
    Aln = []
    for line in Out.split("\n"):
        if line.startswith(">"):
            Aln.append("")
        elif Aln:
            Aln[-1] += line.strip().upper()
    if len(Aln) < 3:
        return None
    A = np.array([list(s) for s in Aln])
    IsBase = np.isin(A, list("ACGT"))
    Ident = []
    for i in range(len(A) - 1):
        Shared = IsBase[i][None, :] & IsBase[i + 1:]
        N = Shared.sum(axis=1)
        Same = ((A[i][None, :] == A[i + 1:]) & Shared).sum(axis=1)
        Ident += list((Same / np.maximum(N, 1))[N >= 30])
    return float(np.median(Ident)) if Ident else None

def Extract(Genome, Hits):
    Regions = ["%s:%d-%d" % (c, s + 1, e) for c, s, e, st, hf, ht in Hits]
    Out = subprocess.run(["samtools", "faidx", Genome] + Regions, capture_output=True, text=True).stdout
    Seqs = []
    for line in Out.split("\n"):
        if line.startswith(">"):
            Seqs.append("")
        elif Seqs:
            Seqs[-1] += line.strip().upper()
    return Seqs

def main():
    try:
        options, remainder = getopt.getopt(sys.argv[1:], 't:L:g:a:E:w:f:d:s:c:h',
            ['tbl=','length=','genome=','anchor=','maxE=','window=','depthFold=','identityGap=','sample=','cpu=','help'])
    except getopt.GetoptError as err:
        print(str(err))
        usage()
        sys.exit(2)

    Anchor = None
    MaxE = 1e-5
    W = 50
    DepthFold = 3.0
    IdGap = 0.15
    Sample = 40
    CPU = 1
    for opt, arg in options:
        if opt in ('-h', '--help'):
            usage()
        elif opt in ('-t', '--tbl'):
            Tbl = arg
        elif opt in ('-L', '--length'):
            L = int(arg)
        elif opt in ('-g', '--genome'):
            Genome = arg
        elif opt in ('-a', '--anchor'):
            Anchor = int(arg)
        elif opt in ('-E', '--maxE'):
            MaxE = float(arg)
        elif opt in ('-w', '--window'):
            W = int(arg)
        elif opt in ('-f', '--depthFold'):
            DepthFold = float(arg)
        elif opt in ('-d', '--identityGap'):
            IdGap = float(arg)
        elif opt in ('-s', '--sample'):
            Sample = int(arg)
        elif opt in ('-c', '--cpu'):
            CPU = int(arg)
    if Anchor is None:
        Anchor = L // 2
    random.seed(1992)

    Hits = ReadHits(Tbl, MaxE)
    Depth = np.zeros(L + 2, dtype=int)
    for c, s, e, st, hf, ht in Hits:
        Depth[hf:ht + 1] += 1
    Starts = list(range(1, L + 1, W))
    WinDepth = [int(np.median(Depth[s:min(s + W, L + 1)])) for s in Starts]
    BaseDepth = max(1, int(np.median([d for d in WinDepth if d > 0])))

    #Identity of windows: the abundant ones, and a baseline from windows at typical depth
    def WindowIdentity(s):
        Cover = [h for h in Hits if h[4] <= s + W // 2 <= h[5]]
        if len(Cover) < 5:
            return None
        Pick = random.sample(Cover, min(Sample, len(Cover)))
        return MedianIdentity(Extract(Genome, Pick), CPU)
    Typical = [s for s, d in zip(Starts, WinDepth) if 0.5 * BaseDepth <= d <= 2 * BaseDepth]
    BaseIds = [i for i in (WindowIdentity(s) for s in random.sample(Typical, min(8, len(Typical)))) if i is not None]
    BaseId = float(np.median(BaseIds)) if BaseIds else 0.5
    print("BASELINE\tdepth %d\tidentity %.3f (%d windows)" % (BaseDepth, BaseId, len(BaseIds)))

    Foreign = []
    for s, d in zip(Starts, WinDepth):
        if d >= DepthFold * BaseDepth:
            i = WindowIdentity(s)
            Young = i is not None and i >= BaseId + IdGap
            print("WINDOW\t%d-%d\tdepth %d (%.1fx)\tidentity %s\t%s" % (s, min(s + W - 1, L), d, d / BaseDepth,
                  "NA" if i is None else "%.3f" % i, "FOREIGN" if Young else "abundant, not young"))
            if Young:
                Foreign.append(s)

    Keep = [1, L]
    Left = [s for s in Foreign if s + W - 1 < Anchor]
    Right = [s for s in Foreign if s > Anchor]
    #Cut at the lowest depth between the anchor and the nearest foreign block on each side
    if Right:
        b = min(Right)
        Gap = range(Anchor, b)
        Keep[1] = min(Gap, key=lambda p: (Depth[p], -p)) if len(Gap) else b - 1
    if Left:
        b = max(Left) + W - 1
        Gap = range(b + 1, Anchor + 1)
        Keep[0] = min(Gap, key=lambda p: (Depth[p], p)) if len(Gap) else b + 1
    print("KEEP\t%d\t%d" % (Keep[0], Keep[1]))

if __name__ == "__main__":
    main()
