#!/usr/bin/env python3

import getopt
import sys
import polars as pl
from Bio import AlignIO

def usage():
        print("Script get the coordinates where alignments end\n")
        print("Usage: CoordinateExtraction.py -i <Alignment> -o <OutputPrefix> [Options]")
        print("\nArguments:")
        print("MANDATORY:")
        print("--input       | -i\t Alignment in fasta format")
        print("--output      | -o\t Output file")
        print("OPTIONAL:")
        print("--BasesOverl  | -B\t Bases to consider side to be aligned (Default: 10)")
        print("--Side        | -S\t Side to consider, can be 'Left', 'Right' or 'Central' (Default: 'Central')")
        print("--help        | -h\t Show this help message and exit")
        exit()

def logmsg(input, output, BasesOverl, Side):
        print("_________________________________________")
        print("Running with the following arguments:")
        print("Input       : ", input)
        print("Output      : ", output)
        print("BasesOverl  : ", BasesOverl)
        print("Side        : ", Side)
        print("_________________________________________")
        print("")

def main():
    try:
        options, remainder = getopt.getopt(sys.argv[1:],'i:o:B:S:h', ['input=','output=','BasesOverl=','Side=','help'])
    except getopt.GetoptError as err:
        print(str(err))
        usage()
        sys.exit(2)

    #Set default values
    BasesOverl = 10
    Side = 'Central'

    #Parse options
    for opt, arg in options:
        if opt in ('-h', '--help'):
            usage()
        elif opt in ('-i', '--input'):
            Input = arg
        elif opt in ('-o', '--output'):
            Output = arg
        elif opt in ('-B', '--BasesOverl'):
            BasesOverl = int(arg)
        elif opt in ('-S', '--Side'):
            Side = arg

    logmsg(Input, Output, BasesOverl, Side)

    #Main code
    ##Read the alignment
    Aln = AlignIO.read( Input , "fasta" )

    ##Create dictionary
    Aln_data = {record.id: list(str(record.seq)) for record in Aln}

    ##Convert the dictionary into a Polars DataFrame
    Aln_DF = pl.DataFrame(Aln_data)
    
    ##Transpose the DF
    Aln_DF = Aln_DF.transpose(include_header=True)

    ##Output the dataframe
    OutputFile = open(Output, 'w')

    #Search for left 
    if Side == 'Left' or Side == 'Central':
        LeftCounter = 1
        NonGapCount = 0
        for i in Aln_DF[0,1:] != '-':
            if i[0] == True:
                NonGapCount += 1
            else:
                NonGapCount = 0
            if NonGapCount >= BasesOverl:
                LeftCounter -= BasesOverl - 1 
                break
            print(LeftCounter)
            LeftCounter+= 1
        #Count only nongap columns
        LeftCounter = (Aln_DF[1,1:LeftCounter+1] != "-").transpose().sum()[0,0]
        #Just a check
        if LeftCounter == 0:
            LeftCounter = 1
        #Output the left side
        OutputFile.write("1\t"+str(LeftCounter)+"\n")


    #Search for right
    if Side == 'Right' or Side == 'Central':
        RightCounter = Aln_DF.shape[1]-1
        NonGapCount = 0
        for i in Aln_DF[0,Aln_DF.columns[::-1][:-1]] != '-':
            if i[0] == True:
                NonGapCount += 1
            else:
                NonGapCount = 0
            if NonGapCount >= BasesOverl:
                RightCounter += BasesOverl -1
                break
            RightCounter-= 1
        #Count only nongap columns
        RightCounter = (Aln_DF[1,1:RightCounter+1] != "-").transpose().sum()[0,0]
        #Output the right side
        OutputFile.write(str(RightCounter)+"\t"+ str((Aln_DF[1,1:] != "-").transpose().sum()[0,0]) +"\n")

    #Close the output file
    OutputFile.close()

if __name__ == "__main__":
    main()
