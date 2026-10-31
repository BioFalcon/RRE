#!/usr/bin/env python3

import getopt
import sys

def usage():
        print("Script to remove the outermost consensus positions of the merged model on one side, in place:")
        print("the alignment columns up to the N-th called column (from the .Columns.txt calls), the matching calls,")
        print("and rows left without bases. The consensus is rewritten from the remaining calls.\n")
        print("Usage: TrimEdge.py -p <Prefix> -s <Side> -n <Positions>")
        print("--prefix      | -p\t Prefix of <Prefix>.aln.fa, <Prefix>.Columns.txt and <Prefix>.Consensus.fa")
        print("--side        | -s\t 'Left' or 'Right'")
        print("--positions   | -n\t Number of consensus positions to remove")
        exit()

def main():
    try:
        options, remainder = getopt.getopt(sys.argv[1:], 'p:s:n:h', ['prefix=','side=','positions=','help'])
    except getopt.GetoptError as err:
        print(str(err))
        usage()
    for opt, arg in options:
        if opt in ('-h', '--help'):
            usage()
        elif opt in ('-p', '--prefix'):
            Prefix = arg
        elif opt in ('-s', '--side'):
            Side = arg
        elif opt in ('-n', '--positions'):
            N = int(arg)

    Calls=open(Prefix+".Columns.txt").read().strip()
    Names=[]; Seqs={}
    for l in open(Prefix+".aln.fa"):
        l=l.strip()
        if l.startswith(">"):
            Names.append(l[1:]); Seqs[l[1:]]=[]
        elif l:
            Seqs[Names[-1]].append(l)
    Seqs={k:"".join(v) for k,v in Seqs.items()}
    ConsName=open(Prefix+".Consensus.fa").readline().strip()

    Called=[i for i,c in enumerate(Calls) if c!="-"]
    if N<=0 or N>=len(Called):
        print("Nothing trimmed (%d positions requested, %d in the model)" % (N, len(Called)))
        return
    if Side=="Left":
        Keep=slice(Called[N-1]+1, len(Calls))
    else:
        Keep=slice(0, Called[-N])
    Calls=Calls[Keep]
    with open(Prefix+".aln.fa","w") as F:
        for n in Names:
            s=Seqs[n][Keep]
            if s.strip("-"):
                F.write(">%s\n%s\n" % (n, s))
    with open(Prefix+".Columns.txt","w") as F:
        F.write(Calls+"\n")
    with open(Prefix+".Consensus.fa","w") as F:
        F.write("%s\n%s\n" % (ConsName, Calls.replace("-","")))
    print("TRIMMED_POSITIONS %d" % N)

if __name__ == "__main__":
    main()
