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
        print("Script to clean aligment and generate a consensus sequence\n")
        print("Usage: MSACuration.py -i <Alignment> -o <OutputPrefix> [Options]")
        print("\nArguments:")
        print("MANDATORY:")
        print("--input        | -i\t Alignment in fasta format")
        print("--output       | -o\t Output prefix")
        print("OPTIONAL:")
        print("--MinCoverage  | -m\t Minimum number of bases for a base to be considered (Default: 5)")
        print("--MainWindow   | -W\t Main Window size used to calculate the threshold used (Default: 100)")
        print("--StepWindow   | -S\t Comma delimited values of the values to use in the stepwise  (Default: 50,25,10,5,2)")
        print("--MinThresh    | -t\t Minimum threshold used to consider a peak as noise (Default: 0.5)")
        print("--SeqID        | -I\t ID of the consensus sequence (Default: Consensus)")
        print("--Plots        | -P\t Generate plots (Default: False)")
        print("--NoiseThresh  | -T\t Threshold (peak value) to consider that alignment has no noise (Default: 0.2)")
        print("--short        | -s\t Short run, enables Plots if enabled (Default: False)")
        print("--zero         | -Z\t Adds a threshold of 0 to the list (Default: False)")
        print("--zeroOnly     | -O\t Only run with threshold zero (Default: False)")
        print("--firstPeak    | -F\t Only remove the first peak, ignoring --NoiseThresh (Default: False)")
        print("--help         | -h\t This beautiful help message :)")
        exit()

def logmsg(input, output, MinCoverage, MainWindow, StepWindow, ID_Sequence, Plots, NoiseThresh, Short, Zero, ZeroOnly, FirstPeak):
        print("_________________________________________")
        print("Running with the following arguments:")
        print("Input       : ", input)
        print("Output      : ", output)
        print("MinCoverage : ", MinCoverage)
        print("MainWindow  : ", MainWindow)
        print("StepWindow  : ", StepWindow)
        print("SeqID       : ", ID_Sequence)
        print("Plots       : ", Plots)
        print("NoiseThresh : ", NoiseThresh)
        print("Short       : ", Short)
        print("Zero        : ", Zero)
        print("ZeroOnly    : ", ZeroOnly)
        print("FirstPeak   : ", FirstPeak)
        print("_________________________________________")
        print("")
            
##Main function
def main():
    try:
        options, remainder = getopt.getopt(sys.argv[1:],'i:o:m:W:S:t:I:PN:ZOFh', ['input=','output=','MinCoverage=','MainWindow=','StepWindow=','MinThresh=','SeqID=','Plots','NoiseThresh=','short','zero','zeroOnly','firstPeak','help'])
    except getopt.GetoptError as err:
        print(str(err))
        usage()
        sys.exit(2)
    
    #Set default values
    ID_Sequence  = "Consensus"
    MainWindow   = 100
    MinCoverage  = 5
    StepWindow   = [25,10,5,2]
    NoiseThresh  = 0.2
    Plots        = False
    Short        = False
    Zero         = False
    ZeroOnly     = False
    FirstPeak    = False

    #Parse arguments
    for opt, arg in options:
        if opt in ('--input','-i'):
            input = arg
        elif opt in ('--output','-o'):
            output = arg
        elif opt in ('--MinCoverage','-m'):
            MinCoverage = int(arg)
        elif opt in ('--MainWindow','-W'):
            MainWindow = int(arg)
        elif opt in ('--StepWindow','-S'):
            StepWindow = [int(i) for i in arg.split(",")]
        elif opt in ('--SeqID','-I'):
            ID_Sequence = arg
        elif opt in ('--Plots','-P'):
            Plots = True
        elif opt in ('--NoiseThresh','-T'):
            NoiseThresh = float(arg)
        elif opt in ('--short','-s'):
            Short = True
            Plots = True
        elif opt in ('--zero','-Z'):
            Zero = True
        elif opt in ('--zeroOnly','-O'):
            ZeroOnly = True
        elif opt in ('--firstPeak','-F'):
            FirstPeak = True
        elif opt in ('--help','-h'):
            usage()

    logmsg(input, output, MinCoverage, MainWindow, StepWindow, ID_Sequence, Plots, NoiseThresh, Short,Zero, ZeroOnly, FirstPeak)

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
    
    ##########Calculate Gaps and Bits of the raw alignment
    print("Step 2: Removimg gappy positions")
    Bits, Gaps = BitGap_Calc(Aln_DF)

    ##########First Plot
    #First Plot before gap removal
    Bits_Sliding = []
    Gaps_Sliding = []
    PlotsLoop=[]
    for i in range(0,len(Bits) - MainWindow):
        Bits_Sliding.append(np.median(Bits[i:i + MainWindow]))
        Gaps_Sliding.append(np.median(Gaps[i:i + MainWindow]))
        #
    if Plots:
        plt.style.use('dark_background')
        Plot1 = plt.figure()
        plt.hist(Bits_Sliding[1:], bins=50)
        plt.xlabel("Bits")
        plt.ylabel("Frequency")
        plt.title("Plot 1 - Bits distribution raw alignment (Sliding Window "+str(MainWindow)+")" )

        Plot2 = plt.figure()
        plt.hist(Gaps_Sliding[1:], bins=50)
        plt.xlabel("Gaps")
        plt.ylabel("Frequency")
        plt.title("Plot 2 - Gaps distribution raw alignment (Sliding Window "+str(MainWindow)+")" )

        Plot3 = plt.figure(figsize=(10, 3))
        plt.subplots_adjust(bottom=0.20)
        plt.plot(range(0,len(Bits_Sliding)), Bits_Sliding, color="#4287f5", lw=0.1)
        plt.xlabel("Position", fontsize=15)
        plt.ylabel("Bits", fontsize=15)
        plt.title("Plot 3 - Bits across raw alignment (Sliding Window "+str(MainWindow)+")", fontsize=18)

        Plot4 = plt.figure(figsize=(10, 3))
        plt.subplots_adjust(bottom=0.20)
        plt.plot(range(0,len(Gaps_Sliding)), Gaps_Sliding, color="#f53b3b", lw=0.1)
        plt.xlabel("Position", fontsize=15)
        plt.ylabel("Gaps", fontsize=15)
        plt.title("Plot 3 - Gaps across raw alignment (Sliding Window "+str(MainWindow)+")", fontsize=18)

    ##########Calculate the threshold
    ##Check if there is at least two different values
    if len(np.unique(Bits)) > 1:
        ##Calculate KDE
        Density=gaussian_kde(Bits, bw_method=0.10)
        DensiVal=Density([x/100 for x in range(-50,100)])

        ##Calculate peaks
        Peaks,_ =find_peaks( DensiVal, prominence=(0.1, None) )
        peakVals = [ [x/100 for x in range(-50,100)][x] for x in Peaks ]

        ##Plot peaks
        if Plots:
            Plot9 = plt.figure()
            plt.plot([x/100 for x in range(-50,100)], DensiVal, color="#4287f5")
            for x in Peaks:
                plt.plot([x/100 for x in range(-50,100)][x], DensiVal[x], "x", color="#f53b3b")
            plt.axvline(x=NoiseThresh, color="#32a852", linestyle="--", label="Threshold")
            plt.xlabel("Bits")
            plt.ylabel("Density")
            plt.title("Plot 9 - Bits density distribution of Bits (Sliding Window "+str(MainWindow)+")")
    
        ##Kill if short#
        if Short:
            pdf = PdfPages(output + ".pdf")
            pdf.savefig(Plot1)
            pdf.savefig(Plot2)
            pdf.savefig(Plot3)
            pdf.savefig(Plot4)
            # pdf.savefig(Plot5)
            # pdf.savefig(Plot6)
            # pdf.savefig(Plot7)
            # pdf.savefig(Plot8)
            pdf.savefig(Plot9)
            pdf.close()
            exit()

        #Add the intact alignment
        if Zero:
            peakVals=[0]+peakVals
            Peaks=np.append(0,Peaks)
        #Zero only overrides others
        if ZeroOnly:
            peakVals=[0]
            Peaks=[0]
        if FirstPeak:
            peakVals=[peakVals[0]]
            Peaks=[Peaks[0]]
            NoiseThresh=peakVals[0]+0.01 #Just for the plot :)
        if any([x <= NoiseThresh for x in peakVals]):
            SelectedPeaksVals = [x for x in peakVals if x <= NoiseThresh]
            SelectedPeaks = [Peaks[peakVals.index(x)] for x in SelectedPeaksVals]
            Slope=[]
            for i in range(1, len(DensiVal)):
                Slope.append( (DensiVal[i]-DensiVal[i-1])/((i/100)-((i-1)/100)) )
            for idx, CurrPeak in enumerate(SelectedPeaks):
                Aln_DFTemp = Aln_DF
                #Calculate Bits
                Bits, Gaps = BitGap_Calc(Aln_DFTemp)
                #Calculate sliding bits
                Bits_Sliding = []
                for i in range(0,len(Bits) - MainWindow):
                    Bits_Sliding.append(np.median(Bits[i:i + MainWindow]))
                if Plots:
                    PlotsLoop.append( plt.figure() )
                    PlotsLoop[idx], axs = plt.subplots(nrows=len(StepWindow)+1, ncols=1, figsize=(16, 8), sharey=True,sharex=True)
                    axs[0].plot(range(0,len(Bits_Sliding)), Bits_Sliding, color="#4a74f0")
                for posSlope in range( SelectedPeaks[idx] , 100):
                    if Slope[posSlope] > 0:
                        break
                posSlopePos = [x/100 for x in range(-50,100)][posSlope]
                for idx2, CurrWindow in enumerate(StepWindow):
                    ##Calculate where to cut 
                    ToKeep=[]
                    for j in range(0,len(Bits) - CurrWindow):
                        if np.median(Bits[j:j + CurrWindow]) > posSlopePos:
                            ToKeep+=list(range(j, j + CurrWindow))
                    ToKeep = set(ToKeep)
                    ToKeep = np.sort( np.array( list(ToKeep) ) + 1 )
                    ##Filter the bits
                    Aln_DFTemp = Aln_DFTemp[ : , [0] + ToKeep.tolist() ]
                    ##Recalculate bits
                    Bits, Gaps = BitGap_Calc(Aln_DFTemp)
                    ##Recaulculate Windowed bits
                    Bits_Sliding = []
                    for i in range(0,len(Bits) - MainWindow):
                        Bits_Sliding.append(np.median(Bits[i:i + MainWindow]))
                    if Plots:
                        axs[idx2+1].plot(range(0,len(Bits_Sliding)), Bits_Sliding, color="#4a74f0")
                ThresholdGaps = Aln_DFTemp.shape[0] - MinCoverage
                toKeep=( np.where([ x < ThresholdGaps for x in Gaps])[0] )
                Aln_DFTemp = Aln_DFTemp[:, [0] + ( toKeep + 1 ).tolist()  ] 
                ConsensusSeq = ""
                for i in Aln_DFTemp.columns[1:]:
                    Bases= list(set(Aln_DFTemp[:,i]))
                    ##Remove gaps
                    if "-" in Bases:
                        Bases.remove("-")
                    ##Count base frequencies
                    BaseFreq = [ Aln_DFTemp.filter(Aln_DFTemp[:,i] == x).shape[0] for x in Bases ]
                    ##Find out if there is a tie
                    if BaseFreq.count(max(BaseFreq)) > 1:
                        ConsensusSeq += "N"
                    else:
                        ConsensusSeq += Bases[BaseFreq.index(max(BaseFreq))]
                File=open(output + ".Consensus__Peak"+str(idx)+".fa", "w")
                File.write(">"+ID_Sequence + "\n")
                for i in range(0, len(ConsensusSeq), 80):
                    File.write(ConsensusSeq[i:i+80] + "\n")
                File.close()
                DF_to_Fasta(Aln_DFTemp, output + "__Peak"+str(idx)+".aln.fa")
        else:
            #Just remove gaps
            ThresholdGaps = Aln_DF.shape[0] - MinCoverage
            toKeep=( np.where([ x < ThresholdGaps for x in Gaps])[0] )
            Aln_DFTemp = Aln_DF[ : , [0] + ( toKeep + 1 ).tolist()  ]
            ConsensusSeq = ""
            for i in Aln_DFTemp.columns[1:]:
                Bases= list(set(Aln_DFTemp[:,i]))
                ##Remove gaps
                if "-" in Bases:
                    Bases.remove("-")
                ##Count base frequencies
                BaseFreq = [ Aln_DFTemp.filter(Aln_DFTemp[:,i] == x).shape[0] for x in Bases ]
                ##Find out if there is a tie
                if BaseFreq.count(max(BaseFreq)) > 1:
                    ConsensusSeq += "N"
                else:
                    ConsensusSeq += Bases[BaseFreq.index(max(BaseFreq))]
            File=open(output + ".Consensus__Peak0.fa", "w")
            File.write(">"+ID_Sequence + "\n")
            for i in range(0, len(ConsensusSeq), 80):
                File.write(ConsensusSeq[i:i+80] + "\n")
            File.close()
            DF_to_Fasta(Aln_DFTemp, output + "__Peak0.aln.fa")
    else:
        #Just remove gaps
        ThresholdGaps = Aln_DF.shape[0] - MinCoverage
        toKeep=( np.where([ x < ThresholdGaps for x in Gaps])[0] )
        Aln_DFTemp = Aln_DF[ : , [0] + ( toKeep + 1 ).tolist()  ]
        ConsensusSeq = ""
        #Obtain consensus sequence
        for i in Aln_DFTemp.columns[1:]:
            Bases= list(set(Aln_DFTemp[:,i]))
            ##Remove gaps
            if "-" in Bases:
                Bases.remove("-")
            ##Count base frequencies
            BaseFreq = [ Aln_DFTemp.filter(Aln_DFTemp[:,i] == x).shape[0] for x in Bases ]
            ##Find out if there is a tie
            if BaseFreq.count(max(BaseFreq)) > 1:
                ConsensusSeq += "N"
            else:
                ConsensusSeq += Bases[BaseFreq.index(max(BaseFreq))]
        File=open(output + ".Consensus__Peak0.fa", "w")
        File.write(">"+ID_Sequence + "\n")
        for i in range(0, len(ConsensusSeq), 80):
            File.write(ConsensusSeq[i:i+80] + "\n")
        File.close()
        DF_to_Fasta(Aln_DFTemp, output + "__Peak0.aln.fa")



    print("Step 5: Generating other outputs")
    #Write outputs
    if Plots:
        pdf = PdfPages(output + ".pdf")
        pdf.savefig(Plot1)
        pdf.savefig(Plot2)
        pdf.savefig(Plot3)
        pdf.savefig(Plot4)
        if "Plot9" in locals():
            pdf.savefig(Plot9)
        if len(PlotsLoop) > 0:
            for i in range(0, len(PlotsLoop)):
                pdf.savefig(PlotsLoop[i])
        pdf.close()

if __name__ == "__main__":
    main()

