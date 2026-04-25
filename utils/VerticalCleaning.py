#!/env/python3

from Bio import AlignIO
import getopt
import sys
import polars as pl
import numpy as np

import matplotlib.pyplot as plt

def usage():
        print("Script to clean Round aligments \n")
        print("Usage: VerticalCleaning.py -i <Alignment> -o <OutputPrefix>")
        print("\nArguments:")
        print("MANDATORY:")
        print("--input        | -i\t Alignment in fasta format")
        print("--output       | -o\t Output prefix")
        print("OPTIONAL:")
        print("--perc         | -P\t Fraction to consider for hand mask (default: 0.1)")
        print("--help         | -h\t This beautiful help message :)")
        exit()

def logmsg(input, output, perc):
        print("_________________________________________")
        print("Running with the following arguments:")
        print("Input       : ", input)
        print("Output      : ", output)
        print("Percentile  : ", perc)
        print("_________________________________________")
        print("")

def main():
    try:
        options, remainder = getopt.getopt(sys.argv[1:],'i:o:P:h', ['input=','output=','perc=','help'])
    except getopt.GetoptError as err:
        print(str(err))
        usage()
        sys.exit(2)

    #Default arguments
    perc = 0.1

    #Parse arguments
    for opt, arg in options:
        if opt in ('--input','-i'):
            input = arg
        elif opt in ('--output','-o'):
            output = arg
        elif opt in ('--perc','-P'):
            perc = float(arg)
        elif opt in ('--help','-h'):
            usage()
            exit(0)

    logmsg(input, output, perc)

    #Main Code
    ##Read alignment into DF    
    Aln = AlignIO.read( input , "fasta" )
    
    ##Create dictionary
    Aln_data = {record.id: list(str(record.seq)) for record in Aln}
    
    ##Convert the dictionary into a Polars DataFrame
    Aln_DF = pl.DataFrame(Aln_data)
    
    ##Transpose the DF
    Aln_DF = Aln_DF.transpose(include_header=True)

    ##Get unique round names
    Rounds=np.unique(np.array([ x.split("___") for x in Aln_DF[Aln_DF.columns[0]] ])[:,1])

    ##Get all names for all sequences
    RoundsAll=np.array([ x.split("___") for x in Aln_DF[Aln_DF.columns[0]] ])[:,1]

    ##Calculate weights for each round
    Weights=np.zeros([Rounds.shape[0],Aln_DF.shape[1]-1])

    for seq in range(0,Aln_DF.shape[0]):
        Name=Aln_DF[seq,0].split("___")[1]
        Weights[np.where(Rounds == Name),:] += 1*(Aln_DF[seq,1:] != "-")


    #Calculate cutoffs based on Percentage of sequences  
    cutoff_array=[]
    for CurrRound in list(range(0,Rounds.shape[0]))[::-1]:
        cutoff_array.append( np.floor(sum(RoundsAll == Rounds[CurrRound])*perc) )
    
    Mask = np.zeros(Aln_DF.shape[1]-1)

    #Create mask based on cutoffs
    Aln_DF = Aln_DF.to_numpy()
    for CurrRound in list(range(0,Rounds.shape[0])):
        RoundCoverage=(Aln_DF[RoundsAll == (Rounds[::-1][CurrRound]),1:] != "-").sum(axis=0)
        Mask[RoundCoverage >= cutoff_array[CurrRound]]+=1

    #Filter alignment
    Aln_DF2 = Aln_DF[:,[0]+list(np.where(Mask >= 1)[0]+1)]

    ##Output alignment in fasta format
    with open(output + ".aln.fa","w") as outfh:
        for row in range(0,Aln_DF2.shape[0]):
            outfh.write(">"+Aln_DF2[row,0]+"\n")
            outfh.write("".join(Aln_DF2[row,1:])+"\n")

if __name__ == "__main__":
    main()