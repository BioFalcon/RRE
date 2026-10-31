#!/usr/bin/env python3

import getopt
import sys

def usage():
        print("Script to test whether a round's new bases are part of the repeat: genomic copies of the model's edge")
        print("(the anchor) must continue colinearly into them. A real end shows as a drop in the fraction of anchored")
        print("copies that continue; drift into a neighbouring repeat or unique flank is supported by few anchored copies.\n")
        print("Usage:")
        print("  AnchoredSupport.py build -r <RoundConsensus> -m <MergedConsensus> -s <Side> -e <Extension> -k <Anchor> -o <Query.fa>")
        print("  AnchoredSupport.py score -q <Query.fa> -t <nhmmer.tblout> -s <Side> -e <Extension> [-a <Span>] [-n <MinCopies>] [-f <MinRatio>]")
        print("\nbuild: query = new bases of the round (first/last Extension bases of its consensus) joined to the outer")
        print("       Anchor bases of the merged model, in model orientation")
        print("score: reads an nhmmer search of the query against the genome, chains colinear hits per locus into copies,")
        print("       counts copies spanning the anchor and how many of them reach each distance into the new bases")
        print("--span        | -a\t Anchored copies must span this many anchor bases next to the new bases (Default: 300)")
        print("--minCopies   | -n\t Anchored copies that must reach a base (Default: 5)")
        print("--minRatio    | -f\t Fraction of the anchored copies that must reach a base (Default: 0.2)")
        print("\nscore prints ANCHORED <copies>, CONTINUING <copies reaching the outermost new base> and")
        print("SUPPORTED <new bases reached by enough anchored copies, counted from the model edge>")
        exit()

def ReadSeq(File):
    return "".join(l.strip() for l in open(File) if not l.startswith(">")).replace("-","").upper()

def Chain(Tbl):
    #Hits on the same target and strand are one copy when they are colinear (small gaps, similar on both sides)
    Hits={}
    for l in open(Tbl):
        if l.startswith("#") or not l.strip():
            continue
        f=l.split()
        Hits.setdefault((f[0],f[11]),[]).append((int(f[4]),int(f[5]),int(f[6]),int(f[7])))
    Copies=[]
    for (Target,Strand),Hs in Hits.items():
        Hs.sort()
        Chains=[]
        for qf,qt,tf,tt in Hs:
            for c in Chains:
                dq=qf-c[1]; dt=(tf-c[3]) if Strand=="+" else (c[3]-tf)
                if -150 < dq < 150 and -150 < dt < 150+0.2*max(dq,0) and abs(dq-dt) < 150:
                    c[1]=max(c[1],qt); c[3]=tt
                    break
            else:
                Chains.append([qf,qt,tf,tt])
        Copies+=[(c[0],c[1]) for c in Chains]
    return Copies

def main():
    if len(sys.argv) < 2 or sys.argv[1] not in ("build","score"):
        usage()
    Mode=sys.argv[1]
    try:
        options, remainder = getopt.getopt(sys.argv[2:], 'r:m:s:e:k:o:q:t:a:n:f:h', ['round=','merged=','side=','extension=',
            'anchor=','output=','query=','tbl=','span=','minCopies=','minRatio=','help'])
    except getopt.GetoptError as err:
        print(str(err))
        usage()

    MinCopies=5
    MinRatio=0.2
    Span=300
    Tolerance=25
    for opt, arg in options:
        if opt in ('-h', '--help'):
            usage()
        elif opt in ('-r', '--round'):
            Round = arg
        elif opt in ('-m', '--merged'):
            Merged = arg
        elif opt in ('-s', '--side'):
            Side = arg
        elif opt in ('-e', '--extension'):
            Extension = int(arg)
        elif opt in ('-k', '--anchor'):
            Anchor = int(arg)
        elif opt in ('-o', '--output'):
            Output = arg
        elif opt in ('-q', '--query'):
            Query = arg
        elif opt in ('-t', '--tbl'):
            Tbl = arg
        elif opt in ('-a', '--span'):
            Span = int(arg)
        elif opt in ('-n', '--minCopies'):
            MinCopies = int(arg)
        elif opt in ('-f', '--minRatio'):
            MinRatio = float(arg)

    if Mode=="build":
        R=ReadSeq(Round); M=ReadSeq(Merged)
        if Side=="Left":
            Seq=R[:Extension]+M[:Anchor]
        else:
            Seq=M[-Anchor:]+R[-Extension:]
        with open(Output,"w") as F:
            F.write(">AnchoredQuery\n%s\n" % Seq)
        return

    QLen=len(ReadSeq(Query)); Anchor=QLen-Extension; Span=min(Span,Anchor)
    Copies=Chain(Tbl)
    #Query layout: Left = [new bases][anchor], Right = [anchor][new bases]; anchored copies span the Span anchor
    #bases next to the new bases (the full anchor is searched so the hits are long enough to be found)
    if Side=="Left":
        Anchored=[a for a,b in Copies if a<=Extension+1+Tolerance and b>=Extension+Span-Tolerance]
        Reach=lambda d: sum(1 for a in Anchored if a<=Extension+1-d+Tolerance)
    else:
        Anchored=[b for a,b in Copies if a<=Anchor-Span+1+Tolerance and b>=Anchor-Tolerance]
        Reach=lambda d: sum(1 for b in Anchored if b>=Anchor+d-Tolerance)
    Need=max(MinCopies, MinRatio*len(Anchored))
    Supported=0
    for d in range(1,Extension+1):
        if Reach(d) < Need:
            break
        Supported=d
    print("ANCHOR_LENGTH %d" % Anchor)
    print("ANCHORED %d" % len(Anchored))
    print("CONTINUING %d" % Reach(Extension))
    print("NEEDED %.1f" % Need)
    print("SUPPORTED %d" % Supported)

if __name__ == "__main__":
    main()
