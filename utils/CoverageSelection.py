#!/usr/bin/env python3

import pandas as pd
import numpy as np
import getopt
import sys
import pybedtools
import pysam
import matplotlib.pyplot as plt
import matplotlib
import statistics

from tqdm import tqdm
pd.options.mode.chained_assignment = None

def usage():
    print("Script to obtain sequences trying to mantain a representative coverage of a consensus\n")
    print("Usage: CoverageSelection.py -i <BlastHits.bed> -o <Prefix> -S <RepeatSize> -T <TargetSeqs>")
    print("\nArguments:")
    print("MANDATORY:")
    print("--input           |-i\t Alignment in fasta format")
    print("--output          |-o\t Output prefix")
    print("--size            |-S\t Consensus size")
    print("--genome          |-G\t Genome fasta file")
    print("OPTIONAL:")
    print("--target          |-T\t Min coverage to aim for (Default:25)")
    print("--percentile      |-P\t Percentile to cutoff (Default:5)")
    print("--window          |-W\t Window to expand the sequence (Default:10000)")
    print("--help            |-h\t This beautiful help message :)")
    exit()

def main():
    try:
        options, remainder = getopt.getopt( sys.argv[1:],
                                        'i:o:S:G:T:P:W:M:C:M:h', 
                                        ['input=','output=','size=',
                                         'target=','percentile=','genome=',
                                         'window=','merge=','percCut','maxSeqs=', 'help'] )
    except getopt.GetoptError as err:
        print(err)
        usage()
        sys.exit(2)

    #Setup default values
    window=10000
    percentile=5
    targetSeqs=25
    mergeSize=500
    percCut=0.40
    LenThresh=1.25
    maxSeqs=100
    for opt, arg in options:
        if opt in ('--input','-i'):
            input = arg
        elif opt in ('--output','-o'):
            output = arg
        elif opt in ('--size','-S'):
            repSize = int(arg)
        elif opt in ('--target','-T'):
            targetSeqs=int(arg)
        elif opt in ('--percentile','-P'):
            percentile=float(arg)
        elif opt in ('--genome','-G'):
            genome=arg
        elif opt in ('--window','-W'):
            window=int(arg)
        elif opt in ('--merge','-M'):
            mergeSize=int(arg)
        elif opt in ('--percCut','-C'):
            percCut=float(arg)
        elif opt in ('--maxSeqs','-M'):
            maxSeqs=int(arg)
        elif opt in ('--help','-h'):
            usage()

    #Read Genome to samfile
    genomeSam = pysam.FastaFile(genome)

    #Read Bed file into memory
    bedBlast = pd.read_table(input,header=None).sort_values(by=4,ascending=False)
    bedBlast.reset_index(inplace=True)
    bedBlast.drop(columns="index",inplace=True)

    #Generate dataframe with only consensus coordinates
    bedConsCoord = pd.DataFrame([[int(i[0]),int(i[1])] for i in [i.split('-') for i in bedBlast[3] ] ])

    #Drop entries that fall bellow threshold
    dropScore=np.percentile( bedBlast[4], percentile )
    bedConsCoord=bedConsCoord[ bedBlast[4] >= dropScore ]
    bedBlast=bedBlast[bedBlast[4] >= dropScore]

    #Generate array to track coverage
    coverageArr=pd.array([0] * repSize)

    #Create DF where output is going to be stored
    finalCoord=pd.DataFrame()
    finalCoordIndiv=pd.DataFrame()

    #Define uncoverable regions
    uncoverRegions=coverageArr.copy()
    for i in range(0,bedBlast.shape[0]):
        tempDomain=bedConsCoord.iloc[i]
        uncoverRegions[tempDomain[0]:tempDomain[1]]-=1

    #Generate copy for plotting purposes
    overallCovArray=abs(uncoverRegions.copy())

    #Remove negative values
    uncoverRegions+=targetSeqs
    uncoverRegions[uncoverRegions<0]=0

    #Setupt an usage index
    useEntry= pd.array([0] * bedBlast.shape[0] )

    #Set up progress bar
    pbar = tqdm(total = bedBlast.shape[0])

    #Open output fasta file
    fastaFile=open(output+".fa",'w')

    #Initialize variables before loop
    coverageAchieved=False
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    currLine=0
    windowCounter=0

    #Main loop
    while coverageAchieved == False:
        #Check if there are still regions to cover
        if (currLine >= bedConsCoord.shape[0]):
            break

        #Decide if chunk is already covered
        if ( uncoverRegions[ bedConsCoord.iloc[currLine][0]:bedConsCoord.iloc[currLine][1]] >= targetSeqs ).any():
            currLine+=1
            pbar.update(1)
            if ( (uncoverRegions >= 25).all() ) or ( currLine == bedConsCoord.shape[0]-1 ):
                break
        else:
            #Extract current entry
            currentEntry = bedBlast.iloc[currLine]

            #Perform window expansion
            expandedEntry = [ currentEntry[0], currentEntry[1]-window, currentEntry[2]+window, currentEntry[5] ]

            #Intersect with Entries
            windowEntry = bedBlast[ ( bedBlast[0] == expandedEntry[0] ) &
                                    ( bedBlast[1] <= expandedEntry[2] ) &
                                    ( bedBlast[2] >= expandedEntry[1] ) &
                                    ( bedBlast[5] == expandedEntry[3] ) ]
            currLine+=1

            #Skip if an entry was already used
            if ( (useEntry[windowEntry.index] == 1).any() ):
                pbar.update(1)
                continue

            #Convert to bedobject
            windowEntry_bedFormat = pybedtools.BedTool.from_dataframe(windowEntry)

            #Sort bed
            windowEntry_bedFormat=windowEntry_bedFormat.sort()

            #Merge coordinates that are close
            windowEntry_bedFormat.merge(d=mergeSize,s=True,o="collapse",c=[4,5,6])
            windowEntry_bedFormat=windowEntry_bedFormat.to_dataframe()

            #Calculate length of sequence
            if( sum(windowEntry_bedFormat["end"] - windowEntry_bedFormat["start"]) > LenThresh*repSize ):
                #Regress to original coordinates
                windowEntry = pd.DataFrame([currentEntry])
                temp = pd.DataFrame([currentEntry])
                temp.columns = windowEntry_bedFormat.columns
                windowEntry_bedFormat = temp
                
            #Add Window number to bed
            windowEntry["Window"]=windowCounter

            #Update index of usage
            useEntry[windowEntry.index]=1

            #Update Final objects
            finalCoordIndiv=pd.concat([finalCoordIndiv,windowEntry])
            finalCoord=pd.concat([finalCoord,pd.DataFrame([list(windowEntry[0])[0],windowEntry[1].min(),windowEntry[2].max(),windowCounter,".",list(windowEntry[5])[0]]).transpose()])

            #Update window counter
            windowCounter+=1

            #Generate Fasta seqeunce
            ##Make header
            fastaHeader=[str(i) for i in [list(windowEntry[0])[0],windowEntry[1].min(),windowEntry[2].max(),list(windowEntry[5])[0],windowCounter] ]
            fastaFile.write(">"+fastaHeader[0]+":"+fastaHeader[1]+"-"+fastaHeader[2]+"("+fastaHeader[3]+")__Window"+fastaHeader[4]+"\n")
            
            ##Obtain sequence
            sequence=""
            for i in windowEntry_bedFormat.index:
                sequence+=genomeSam.fetch(windowEntry_bedFormat["chrom"][i],int(windowEntry_bedFormat["start"][i]),int(windowEntry_bedFormat["end"][i]))
            
            #Change to uppercase
            sequence=sequence.upper()
            
            ##Reverse complement seq if needed
            if windowEntry_bedFormat.iloc[0]["strand"].split(",")[0] == "-":
                sequence="".join(complement.get(base, base) for base in reversed(sequence))

            ##Write sequence into file
            for i in range(0, len(sequence), 80):
                fastaFile.write(sequence[i:i + 80] + '\n')

            #Update Coverage
            for i in range(0,windowEntry_bedFormat.shape[0]):
                tempRange=[int(i) for i in windowEntry_bedFormat.iloc[i]['name'].split("-")]
                tempRange=range(tempRange[0],tempRange[1])
                coverageArr[tempRange]+=1
                uncoverRegions[tempRange]+=1

            #Update progress bar
            pbar.update(1)

            #Decide if cycle should break
            if ( (uncoverRegions >= targetSeqs).all() ) or ( currLine >= bedConsCoord.shape[0]-1 ) or ( maxSeqs == windowCounter ):
                break

    #Close fasta file
    fastaFile.close()

    #Save value for consensus cleaning 
    CutCoverageFile=open("./"+output+".CutCoverage.txt",'w')
    CutCoverageFile.write( str(round(statistics.median(coverageArr[coverageArr !=0] )*percCut)) )
    CutCoverageFile.close()

    #Write bed files 
    finalCoordIndiv.to_csv("./"+output+".IndivHits.bed", sep='\t', header=False, index=False)
    finalCoord.to_csv("./"+output+".Window.bed", sep='\t', header=False, index=False)

    #Make plots
    #Make plots for coverage
    fig = plt.figure()
    plt.style.use('dark_background')
    plt.subplots_adjust(hspace=0.5)
    ax1 = fig.add_subplot(2,1,2)
    ax1.grid(color='#cccccc', linestyle='--', linewidth=0.5)
    ax1.plot(range(0,len(coverageArr)),coverageArr,color="#4287f5")
    ax1.get_xaxis().set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))
    ax1.get_yaxis().set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))
    ax1.title.set_text('Subset Coverage')
    ax2 = fig.add_subplot(2,1,1)
    ax2.grid(color='#cccccc', linestyle='--', linewidth=0.5)
    ax2.plot(range(0,len(overallCovArray)),overallCovArray,color="#f54242")
    ax2.title.set_text('Overall Coverage')
    ax2.get_xaxis().set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))
    ax2.get_yaxis().set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))
    fig.supylabel("Coverage (seqs)")
    fig.supxlabel("Position (bp)")
    fig.savefig(output+"_CoveragePlot.pdf", format="pdf")


if __name__ == "__main__":
    main()
