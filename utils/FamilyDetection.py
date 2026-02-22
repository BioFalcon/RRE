#!/usr/bin/env python3

from Bio import AlignIO
from scipy.stats import gaussian_kde
from scipy.signal import find_peaks
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from scipy.spatial.distance  import squareform
from matplotlib.backends.backend_pdf import PdfPages
import sys
import getopt
import numpy as np
import matplotlib.pyplot as plt

def usage():
        print("Script to determine if multiple families are present in an alignment\n")
        print("Usage: FamilyDetection.py -i <Alignment> -o <OutputPrefix> [Options]")
        print("\nArguments:")
        print("MANDATORY:")
        print("--input        | -i\t Alignment in fasta format")
        print("OPTIONAL:")
        print("--output       | -o\t Output prefix (Default: Output)")
        print("--Peakdist     | -d\t Distance between peaks to cluster (Default: 0.18)")
        print("--help         | -h\t This beautiful help message :)")
        print("\n")

def logmsg(input, output, Peakdist):
        print("_________________________________________")
        print("Running with the following arguments:")
        print("Input       : ", input)
        print("Output      : ", output)
        print("Peakdist    : ", Peakdist)
        print("_________________________________________")
        print("")

def main():
    try:
        options, remainder = getopt.getopt(sys.argv[1:],'i:o:d:h', ['input=','output=','Peakdist=','help'])
    except getopt.GetoptError as err:
        print(str(err))
        usage()
        sys.exit(2)

    #Set defaults
    Peakdist = 0.18
    output = "Output"
    Plots=True

    #Parse arguments
    for opt, arg in options:
        if opt in ('--input','-i'):
            input = arg
        elif opt in ('--output','-o'):
            output = arg
        elif opt in ('--Peakdist','-d'):
            Peakdist = float(arg)
        elif opt in ('--help','-h'):
            usage()
            sys.exit()

    logmsg(input, output, Peakdist)

    #Main code
    ##Read alignment
    Aln = AlignIO.read(input, "fasta")

    ##Convert to dictionary
    Aln_data = {record.id: list(str(record.seq)) for record in Aln}

    ##Get keys
    Aln_keys = list(Aln_data.keys())

    ##Convert to numpy array
    Aln_DF = np.array([Aln_data[i] for i in Aln_keys])

    ##Calculate pairwise identity
    Identity = np.empty([ Aln_DF.shape[0] , Aln_DF.shape[0] ])
    for i in range(Aln_DF.shape[0]):
        for j in range(Aln_DF.shape[0]):
            Selection = (Aln_DF[i] != "-") & (Aln_DF[j] != "-")
            if any(Selection):
                Identity[i,j] = sum((Aln_DF[i,Selection] == Aln_DF[j,Selection]))/sum(Selection)
            else:
                Identity[i,j] = 0
                
    ##Covert to distance
    Distance = 1 - Identity
    
    #Generate density function
    ##Check if all sequences are the same
    if len(np.unique(Distance[np.tril(Distance) != 0])) > 1:
        Density = gaussian_kde(Distance[np.tril(Distance) != 0])
        DensiVal = Density([x/100 for x in range(-50,150)])

        ##Find peaks
        Peaks,_ = find_peaks( DensiVal, prominence=(0.40, None) )
        peakVals = [ [x/100 for x in range(-50,150)][x] for x in Peaks ]

        #Generate density plot
        Plot1 = plt.figure()
        plt.plot([x/100 for x in range(-50,150)], DensiVal, color="#4287f5", linewidth=3)
        for x in Peaks:
            plt.plot([x/100 for x in range(-50,150)][x], DensiVal[x], "x", color="#f53b3b", markersize=10, markeredgewidth=3)
        plt.xlabel("Distance", fontsize=15)
        plt.ylabel("Density", fontsize=15)
        plt.xlim(0, 1)
        plt.title("Plot 1 - Density of pairwise identity", fontsize=18)

        ##Select peaks if there are clusters
        SelectedPeaks = np.array([])
        if len( Peaks ) > 1:
            #Make all vs all peak comparisons
            for i in range(len(Peaks)):
                for j in range(i+1, len(Peaks)):
                    if abs(peakVals[i]-peakVals[j]) > Peakdist:
                        #Calculate midpoint and store
                        SelectedPeaks = np.append(SelectedPeaks, (peakVals[i]+peakVals[j])/2 )

        ##Remove duplicates
        SelectedPeaks = np.unique(SelectedPeaks)
        
        #Sort by value
        SelectedPeaks = np.sort(SelectedPeaks)

        ##Continue only if there are selected peaks
        if len(SelectedPeaks) > 0:
            ##Convert to condensed distance matrix
            CondenDis = squareform(Distance)

            ##Generate linkage matrix
            LinkageM = linkage(CondenDis, 'average')

            #Add to plot cutoffs
            plt.axvline(x=SelectedPeaks[0], color='r', linestyle='--', label=f"Cutoff = {SelectedPeaks[0]}")

            ##Cluster sequences
            clusters = fcluster(LinkageM, t=SelectedPeaks[0], criterion='distance')

            #Plot dendrogram
            Plot2 = plt.figure(figsize=(10, 5))
            dendrogram(LinkageM, labels=range(1, len(Identity) + 1))
            plt.axhline(y=SelectedPeaks[0], color='r', linestyle='--', label=f"Cutoff = {SelectedPeaks[0]}")
            plt.title("Dendrogram with Cutoff for Clustering")
            plt.xlabel("Sequence Index")
            plt.ylabel("Dissimilarity")
            plt.legend()
        else:
            clusters = np.zeros(Aln_DF.shape[0])
    else:
        Plots=False
        clusters = np.zeros(Aln_DF.shape[0])

    ##Save plots
    if Plots:
        pdf = PdfPages(output + ".Plots.pdf")
        pdf.savefig(Plot1)        
        if len(SelectedPeaks) > 0:
            pdf.savefig(Plot2)
        pdf.close()

    #Output Alignment
    Family = 0
    for i in np.unique(clusters):
        OutAlnFile = open(output + ".Cluster_"+str(Family)+".fa.aln", "w")
        CurrentSeqInCluster = np.array(Aln_keys)[clusters == i]
        for j in CurrentSeqInCluster:
            OutAlnFile.write(">"+j+"\n")
            for i in range(0, len(Aln_data[j]), 80):
                OutAlnFile.write("".join(Aln_data[j][i:i+80]) + "\n")
        OutAlnFile.close()
        Family += 1
        
if __name__ == "__main__":
    main()