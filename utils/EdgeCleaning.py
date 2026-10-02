#!/usr/bin/env python3

from Bio import AlignIO
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import gaussian_kde
from scipy.signal import find_peaks
import getopt
import sys
import polars as pl
import matplotlib.pyplot as plt
import numpy as np

##Functions
def DF_to_Fasta(DF, OutFile):
    File=open(OutFile, "w")
    for i in range(0, DF.shape[0]):
        File.write(">"+ DF.row(i)[0]+ "\n")
        Sequence="".join(DF.row(i)[1:])
        for j in range(0, len(Sequence), 80):
            File.write(Sequence[j:j+80] + "\n")
    File.close()

def BitGap_Calc(DF):
    Gaps=[]
    Bits=[]
    Lmda = np.log2(min(5, DF.shape[0]))
    for i in DF.columns[1:]:
        CurrCol = DF[i].value_counts().sort(by="count")
        YesGaps = CurrCol.filter(CurrCol[:,0] == "-")
        if len(YesGaps) > 0:
            Gaps.append( CurrCol.filter(CurrCol[:,0] == "-")[0,1] )
            GapFrac = CurrCol.filter(CurrCol[:,0] == "-")[0,1]/CurrCol["count"].sum()
        else:
            Gaps.append(0)
            GapFrac = 0
        #Calculate bits
        ##Calculate Probabilities
        NoGaps = CurrCol.filter(CurrCol[:,0] != "-")
        Probs = NoGaps["count"]/NoGaps["count"].sum()
        Tx = -1*(Probs*np.log2(Probs)).sum()*pow(Lmda, -1)
        Ctrident = pow(( 1 - Tx), 2)*pow(( 1 - GapFrac),0.5)
        Bits.append( Ctrident )
        ##Calculate entropy
        #Entropy = -1*(Probs*np.log2(Probs)).sum()
        ##Calculate bits
        #Bits.append(np.log2(5) - Entropy)
    return Bits, Gaps

def usage():
        print("Script to remove overextension from the ends of a repeat\n")
        print("Usage: EdgeCleaning.py -i <Alignment> -o <OutputPrefix> [Options]")
        print("\nArguments:")
        print("MANDATORY:")
        print("--input        | -i\t Alignment in fasta format")
        print("--output       | -o\t Output prefix")
        print("OPTIONAL:")
        print("--winSize      | -w\t (Default: 10)")
        print("--limit        | -l\t (Default: 125)")
        print("--threshold    | -t\t (Default: 0.05)")
        print("--minSeqs      | -m\t (Default: 5)")
        print("--seqID        | -s\t (Default: Sequence)")
        print("--help         | -h\t This beautiful help message :)")
        exit()

def logmsg(input, output, winSize, limit, threshold, minSeqs, seqID):
        print("_________________________________________")
        print("Running with the following arguments:")
        print("Input       : ", input)
        print("Output      : ", output)
        print("Window size : ", winSize)
        print("Limit       : ", limit)
        print("Threshold   : ", threshold)
        print("Minimum Seqs: ", minSeqs)
        print("Sequence ID : ", seqID)
        print("_________________________________________")
        print("")

##Main function
def main():
    try:
        options, remainder = getopt.getopt(sys.argv[1:],'i:o:w:l:t:m:s:h', ['input=','output=','winSize=','limit=','threshold=','minSeqs=','seqID=','help'])
    except getopt.GetoptError as err:
        print(str(err))
        usage()
        sys.exit(2)

    #Set default values
    winSize   = 10
    limit     = 125
    threshold = 0.05
    minSeqs   = 5
    seqID     = "Sequence"
    
    #Parse arguments
    for opt, arg in options:
        if opt in ('--input','-i'):
            input = arg
        elif opt in ('--output','-o'):
            output = arg
        elif opt in ('--winSize','-w'):
            winSize = int(arg)
        elif opt in ('--limit','-l'):
            limit = int(arg)
        elif opt in ('--threshold','-t'):
            threshold = float(arg)
        elif opt in ('--minSeqs','-m'):
            minSeqs = int(arg)
        elif opt in ('--seqID','-s'):
            seqID = arg
        elif opt in ('--help','-h'):
            usage()

    logmsg(input, output, winSize, limit, threshold, minSeqs, seqID)

    #Main Code
    ##########Read alignment into DF
    print("Step 1: Reading alignment")
    Aln = AlignIO.read( input , "fasta" )
    ##Create dictionary
    Aln_data = {record.id: list(str(record.seq)) for record in Aln}
    ##Convert the dictionary into a Polars DataFrame
    Aln_DF = pl.DataFrame(Aln_data)
    ##Transpose the DF
    Aln_DF = Aln_DF.transpose(include_header=True)

    #Remove gappy regions from the alignment
    Bits, Gaps = BitGap_Calc(Aln_DF)
    Aln_DF = Aln_DF[:, [True] + (np.array([Aln_DF.shape[0]-x for x in Gaps]) >= minSeqs).tolist()]

    ##########Calculate Gaps and Bits of the raw alignment
    print("Step 2: Filtering")
    Bits, Gaps = BitGap_Calc(Aln_DF)

    #Calculate sliding window for bits
    SlidingBits = []
    for i in range(0,len(Bits) - winSize + 1):
        SlidingBits.append(np.median(Bits[i:i + winSize]))

    SlidingDerivative = np.diff(SlidingBits)

    #Calculate threshold from left size
    ##Calculate minimum under threshold
    if winSize < limit:
        LeftCut = True
        LeftMinimum = np.argmin(SlidingBits[0:(limit-winSize)])
        if SlidingBits[0:(limit-winSize)][LeftMinimum] > threshold:
            LeftMinimum = 0
            LeftCut = False

        LeftMinSliding=LeftMinimum
        while LeftCut:
            if SlidingBits[LeftMinSliding] < threshold:
                LeftMinSliding += 1
            else:
                if SlidingDerivative[LeftMinSliding-1] <=0:
                    LeftMinSliding += 1
                else:
                    LeftMinSliding = LeftMinSliding + np.ceil(winSize/2)
                    LeftCut = False
            if LeftMinSliding >= (limit-winSize):
                LeftMinSliding = 0
                LeftCut = False

        #Calculate threshold from right size
        RightCut = True
        RightMinimum = np.argmin(SlidingBits[-(limit-winSize):]) + (len(SlidingBits) - (limit-winSize))
        if SlidingBits[RightMinimum] > threshold:
            RightMinimum = len(SlidingBits) + winSize
            RightCut = False

        RightMinSliding=RightMinimum
        while RightCut:
            if SlidingBits[RightMinSliding] < threshold:
                RightMinSliding -= 1
            else:
                if SlidingDerivative[RightMinSliding-1] >=0:
                    RightMinSliding -= 1
                else:
                    RightMinSliding = RightMinSliding +np.ceil(winSize/2)
                    RightCut = False
            if RightMinSliding <= len(SlidingBits) - (limit-winSize):
                RightMinSliding = len(SlidingBits) + winSize
                RightCut = False
    else:
        LeftMinSliding = 0
        RightMinSliding = len(SlidingBits) + winSize
        
    #Create plot with thresholds and cutoffs
    pdf = PdfPages(output + ".pdf")
    Plot= plt.figure()
    plt.plot(range(0,len(SlidingBits)), SlidingBits)
    plt.axhline(y=threshold, color='r', linestyle='--')
    plt.axvline(x=LeftMinSliding, color='g', linestyle='--')
    plt.axvline(x=RightMinSliding, color='g', linestyle='--')
    plt.xlabel("Position (Sliding window)")
    plt.ylabel("Trident Bits")
    plt.title("Trident Bits across the alignment")
    pdf.savefig(Plot)
    pdf.close()

    #Filter alignment
    Filtered_DF = Aln_DF[:, [0] + list(range(int(LeftMinSliding+1), int(RightMinSliding)))]


    #Output filtered alignment
    DF_to_Fasta(Filtered_DF, output + ".aln.fa")

    #Make consensus sequence
    ConsensusSeq = ""
    for i in Filtered_DF.columns[1:]:
        Bases= list(set(Filtered_DF[:,i]))
        ##Remove gaps
        if "-" in Bases:
            Bases.remove("-")
        ##Count base frequencies
        BaseFreq = [ Filtered_DF.filter(Filtered_DF[:,i] == x).shape[0] for x in Bases ]
        ##Find out if there is a tie
        if BaseFreq.count(max(BaseFreq)) > 1:
            ConsensusSeq += "N"
        else:
            ConsensusSeq += Bases[BaseFreq.index(max(BaseFreq))]
    File=open(output + ".Consensus.fa", "w")
    File.write(">"+ seqID + "\n")
    for i in range(0, len(ConsensusSeq), 80):
        File.write(ConsensusSeq[i:i+80] + "\n")
    File.close()

if __name__ == "__main__":
    main()