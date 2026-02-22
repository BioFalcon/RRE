process HEEA_Extension{
    label "HEEA_Extension"

    input:
    file(Consensus)
    file(HMMRDB)
    file(ConsensusFasta)
    file(GenomeIndex)
    file(HMMROut)
    file(outDir)
    file(Genome)
    file(CurateScript)
    file(OutlierScript)
    file(CoordExtractScript)

    output:
    file("FINAL_HEEA.CHECK")
    
    script:
    """
    #Load ID of the repeat into a variable
    Repeat=\$( cat ${Consensus} )
    FullID=\$( grep "\${Repeat}#" ${ConsensusFasta} | sed 's/>//;s/ .*//' )

    #Save nf directory into a variable
    oldWorkDir=\$(pwd)

    #Move to new WorkDir 
    cd ${outDir}/WorkDir
    mkdir -p \${Repeat}/HEEA
    cd \${Repeat}/HEEA

    #Backup new path
    NewPath=\$(pwd)
    #Establish how many sequences will be used
    SampleSeqs=25

    #Initialize variables
    RoundHEEA=01
    GoodCheck="Yes"
    SkipRounds=No
    ExtensionSize=${params.extension}
    Increase=\$(echo "\${ExtensionSize}" | awk '{\$0=\$0*0.25;printf("%d\\n",\$0+=\$0<0?0:0.9)}' )

    #Init checks
    FailedHEEA=No
    FinalFiles=No
    ExtendLeft=True
    ExtendRight=True

    #Check if Staging dir exist
    if [ ! -d Staging ] || [ ! -f ./Staging/GoodCHECK ] ;then
        mkdir -p Staging
        ln -s -f \${oldWorkDir}/* ./Staging/

        #Extract the blast entries for that repeat
        awk -v Repeat="\${Repeat}" '{OFS="\\t"} \$3==Repeat{print}' ./Staging/${HMMROut} > ./Staging/HMMROut_Subset.txt

        #Make check for Staging
        > ./Staging/GoodCHECK
    fi

    #Make the first round
    if [ ! -d Round_00 ] || [ ! -f ./Round_00/GoodCHECK ]; then
        #Create dir
        mkdir -p Round_00
        
        #Select which sequences will be used 
        sort -k14,14nr ./Staging/HMMROut_Subset.txt | \
        cut -f1,9,10,12,14 | \
        awk '{OFS="\\t"}{ if (\$4 == "+"){print \$1,\$2,\$3,".",\$5,\$4}else{print \$1,\$3,\$2,".",\$5,\$4}}' | \
        sort -k5,5nr > ./Round_00/02_CurrentConsensi.Round00.Extended.bed 

        #Extract previous consensus
        samtools faidx ./Staging/${ConsensusFasta} \${FullID} > ./Round_00/07_CurrentConsensi.Round00.Extended.Curated.Consensus.fa

        #Make check
        > ./Round_00/GoodCHECK
    fi

    #Check if other Rounds exist and backup if needed
    if [ \$(ls -d Round_*| wc -l) -gt 1 ];then
        #Check if Final dir is there
        if [ -d Final ] ;then
            #Set variable to not do rounds
            SkipRounds=Yes
            #Check if Final is OK
            if [ -f ./Final/GoodCHECK ];then
                FinalFiles=Yes
            fi
        else
            #Get which is the last round, and back it up
            RoundHEEA=\$(ls -d ./Round_*| tail -n 1| sed 's/.\\/Round_//')

            #Check if it is complete and increase round
            if [ -f ./Round_\${RoundHEEA}/GoodCHECK ];then
                RoundHEEA=\$(echo \${RoundHEEA} | awk '{printf "%02d\\n", \$1+1}')
            else
                #Make backup of failed round 
                BackupNum=01
                Back=True
                while [ \${Back} == "True" ];do
                    if [ -d ./BackUp_Round\${RoundHEEA}_\${BackupNum} ];then
                        BackupNum=\$(echo \${BackupNum} | awk '{printf "%02d\\n", \$1+1}')
                    else
                        mv ./Round_\${RoundHEEA} ./BackUp_Round\${RoundHEEA}_\${BackupNum}
                        Back=False
                    fi
                done
            fi
        fi
    fi

    #Da main loop
    while [ \${SkipRounds} == "No" ]; do
        #Make dir for this round
        mkdir Round_\${RoundHEEA}

        #Set previous round
        PrevRoundHEEA=\$(echo \${RoundHEEA} | awk '{printf "%02d\\n", \$1-1}')

        #Extract coordinates in sense
        ln -s -f \${NewPath}/Round_\${PrevRoundHEEA}/02_CurrentConsensi.Round\${PrevRoundHEEA}.Extended.bed ./Round_\${RoundHEEA}/00_CurrentConsensi.Round\${RoundHEEA}.coord.bed

        #Extend coordinates
        if [ \${ExtendRight} == "True" ] && [ \${ExtendLeft} == "True" ];then
            bedtools slop \\
                -b \${ExtensionSize} \\
                -i ./Round_\${RoundHEEA}/00_CurrentConsensi.Round\${RoundHEEA}.coord.bed \\
                -g <(cut -f1,2 ./Staging/${GenomeIndex}) \\
                > ./Round_\${RoundHEEA}/01_CurrentConsensi.Round\${RoundHEEA}.Extended.bed
        elif [ \${ExtendRight} == "True" ] && [ \${ExtendLeft} == "False" ];then
            bedtools slop \\
                    -r \${ExtensionSize} \\
                    -l 0 \\
                    -i ./Round_\${RoundHEEA}/00_CurrentConsensi.Round\${RoundHEEA}.coord.bed \\
                    -g <(cut -f1,2 ./Staging/${GenomeIndex}) \\
                    > ./Round_\${RoundHEEA}/01_CurrentConsensi.Round\${RoundHEEA}.Extended.bed
        elif [ \${ExtendRight} == "False" ] && [ \${ExtendLeft} == "True" ];then
            bedtools slop \\
                -l \${ExtensionSize} \\
                -r 0 \\
                -i ./Round_\${RoundHEEA}/00_CurrentConsensi.Round\${RoundHEEA}.coord.bed \\
                -g <(cut -f1,2 ./Staging/${GenomeIndex}) \\
                > ./Round_\${RoundHEEA}/01_CurrentConsensi.Round\${RoundHEEA}.Extended.bed
        fi

        #Select the sequences to be used
        CounterSeqs=0
        FileLine=1
        TotalLines=\$(wc -l < ./Round_\${RoundHEEA}/01_CurrentConsensi.Round\${RoundHEEA}.Extended.bed)
        > ./Round_\${RoundHEEA}/02_CurrentConsensi.Round\${RoundHEEA}.Extended.bed

        while [ \$CounterSeqs -lt \${SampleSeqs} ];do
            #Extract current line
            head -n \${FileLine} ./Round_\${RoundHEEA}/01_CurrentConsensi.Round\${RoundHEEA}.Extended.bed | tail -n 1  > ./Round_\${RoundHEEA}/Temp.bed

            #Check if there is no overlap
            if [[ \$(bedtools intersect -a ./Round_\${RoundHEEA}/Temp.bed -b ./Round_\${RoundHEEA}/02_CurrentConsensi.Round\${RoundHEEA}.Extended.bed | wc -l) -eq 0 ]] ;then
                cat ./Round_\${RoundHEEA}/Temp.bed >> ./Round_\${RoundHEEA}/02_CurrentConsensi.Round\${RoundHEEA}.Extended.bed
                CounterSeqs=\$((\${CounterSeqs}+1))
            fi
        
            #Increase FileLine
            FileLine=\$((\${FileLine}+1))

            if [[ \${FileLine} -gt \${TotalLines} ]];then
                break
            fi
        done

        #Check if it has enough sequences to continue
        if [ \$(wc -l < ./Round_\${RoundHEEA}/02_CurrentConsensi.Round\${RoundHEEA}.Extended.bed ) -lt 10 ] ;then
            break
        fi

        #Get fasta sequences
        bedtools getfasta \\
            -fi ./Staging/${Genome} \\
            -fo ./Round_\${RoundHEEA}/03_CurrentConsensi.Round\${RoundHEEA}.Extended.fa \\
            -bed ./Round_\${RoundHEEA}/02_CurrentConsensi.Round\${RoundHEEA}.Extended.bed \\
            -s

        #Make alignment
        mafft \\
            --localpair \\
            --maxiterate 1000 \\
            --thread ${task.cpus} \\
            ./Round_\${RoundHEEA}/03_CurrentConsensi.Round\${RoundHEEA}.Extended.fa \\
            > ./Round_\${RoundHEEA}/04_CurrentConsensi.Round\${RoundHEEA}.Extended.aln \\
            2> ./Round_\${RoundHEEA}/04_CurrentConsensi.Round\${RoundHEEA}.Extended.MAFFTlog

        #Extract and uppercase 
        seqkit seq \\
            -u \\
            ./Round_\${RoundHEEA}/04_CurrentConsensi.Round\${RoundHEEA}.Extended.aln \\
            > 04_TEMP.aln
        mv 04_TEMP.aln ./Round_\${RoundHEEA}/04_CurrentConsensi.Round\${RoundHEEA}.Extended.aln

        #Curate alignment and make consensus
        ./Staging/${CurateScript} \\
            --input ./Round_\${RoundHEEA}/04_CurrentConsensi.Round\${RoundHEEA}.Extended.aln \\
            --output ./Round_\${RoundHEEA}/05_CurrentConsensi.Round\${RoundHEEA}.Extended.Curated \\
            --SeqID \${Repeat}#___Round\${RoundHEEA} \\
            --Plots \\
            --NoiseThresh ${params.noiseThreshold}

        #Select cleanest
        NumberOfPeaks=\$(ls ./Round_\${RoundHEEA}/05_CurrentConsensi.Round\${RoundHEEA}.Extended.Curated__Peak*.aln.fa | wc -l)
        
        SelectedPeak=\$(( \${NumberOfPeaks} - 1 ))
        HighestMedian=0


        echo "Selected Peak: \${SelectedPeak}"
        echo "Highest Median: \${HighestMedian}"

        #Copy alignment with the best score
        cp ./Round_\${RoundHEEA}/05_CurrentConsensi.Round\${RoundHEEA}.Extended.Curated.Consensus__Peak\${SelectedPeak}.fa ./Round_\${RoundHEEA}/07_CurrentConsensi.Round\${RoundHEEA}.Extended.Curated.Consensus.fa
        cp ./Round_\${RoundHEEA}/05_CurrentConsensi.Round\${RoundHEEA}.Extended.Curated__Peak\${SelectedPeak}.aln.fa ./Round_\${RoundHEEA}/07_CurrentConsensi.Round\${RoundHEEA}.Extended.Curated.aln.fa

        #Align the new consensus to the previous one
        mafft \
            --localpair \
            --maxiterate 1000 \
            --thread ${task.cpus} \
            <( cat ./Round_\${PrevRoundHEEA}/07_CurrentConsensi.Round\${PrevRoundHEEA}.Extended.Curated.Consensus.fa \
               ./Round_\${RoundHEEA}/07_CurrentConsensi.Round\${RoundHEEA}.Extended.Curated.Consensus.fa ) \
            > ./Round_\${RoundHEEA}/08_CurrentConsensi.Round\${RoundHEEA}.Consensi.aln.fa
        
        #Extract Coordinates
        ./Staging/${CoordExtractScript} \
            --input ./Round_\${RoundHEEA}/08_CurrentConsensi.Round\${RoundHEEA}.Consensi.aln.fa \
            --output ./Round_\${RoundHEEA}/08_CurrentConsensi.Round\${RoundHEEA}.Consensi.coord.bed \
            --Side Central

        ##Get Length of extension into variables
        LengthExtension=\$(sed '/>/d' ./Round_\${RoundHEEA}/07_CurrentConsensi.Round\${RoundHEEA}.Extended.Curated.Consensus.fa | tr '\\n' ' '| sed 's/ //g'| wc -m)
        ConsensusID=\$(grep ">" ./Round_\${RoundHEEA}/07_CurrentConsensi.Round\${RoundHEEA}.Extended.Curated.Consensus.fa | sed 's/>//')
        LeftSide_Coord=\$( head -n 1 ./Round_\${RoundHEEA}/08_CurrentConsensi.Round\${RoundHEEA}.Consensi.coord.bed | cut -f2 )
        RightSide_Coord=\$( tail -n 1 ./Round_\${RoundHEEA}/08_CurrentConsensi.Round\${RoundHEEA}.Consensi.coord.bed | cut -f1 )
        PreviousLength=\$(sed '/>/d' ./Round_\${PrevRoundHEEA}/07_CurrentConsensi.Round\${PrevRoundHEEA}.Extended.Curated.Consensus.fa | tr '\\n' ' '| sed 's/ //g'| wc -m)

        if [ \${LengthExtension} -gt \${PreviousLength} ];then
            ##Check Left
            if [ \${LeftSide_Coord} -lt \${Increase} ]; then
                ExtendLeft=False
            fi

            ##Check Right
            if [ \$(( \${LengthExtension} - \${RightSide_Coord} )) -lt \${Increase} ]; then
                ExtendRight=False
            fi

        else
            ExtendRight=False
            ExtendLeft=False
        fi

        #If it made it this far, make check 
        > ./Round_\${RoundHEEA}/GoodCHECK

        if [ \${ExtendLeft} == "False" ] && [ \${ExtendRight} == "False" ];then
            break
        fi
    
        #Break if consensus is longer than N bases
        if [ \$( seqkit fx2tab --length --name ./Round_\${RoundHEEA}/07_CurrentConsensi.Round\${RoundHEEA}.Extended.Curated.Consensus.fa | cut -f2)  -gt 25000 ];then
            break
        fi

        #Increase Round number
        RoundHEEA=\$(echo \${RoundHEEA} | awk '{printf "%02d\\n", \$1+1}')
    done

    #Wrap-up
    if [ \${FinalFiles} == "No" ];then
        #Make Final directory
        mkdir -p ./Final

        #Determine last OK round
        LastSuccRound=\$(ls ./Round_*/GoodCHECK| tail -1| sed 's/\\/GoodCHECK//;s/.*Round_//')

        if [ \${LastSuccRound} == '00' ];then
            echo ">\${Repeat}#___Failed" > ./Final/FinalConsensi.aln.fa
            echo "A" >> ./Final/FinalConsensi.aln.fa
            echo ">\${Repeat}#___Failed" > ./Final/FinalConsensi.consensus.fa
            echo "A" >> ./Final/FinalConsensi.consensus.fa
            FailedHEEA=Yes
        else
            cp ./Round_\${LastSuccRound}/07_CurrentConsensi.Round\${LastSuccRound}.Extended.Curated.aln.fa       ./Final/FinalConsensi.aln.fa
            cp ./Round_\${LastSuccRound}/07_CurrentConsensi.Round\${LastSuccRound}.Extended.Curated.Consensus.fa ./Final/FinalConsensi.consensus.fa
        fi

        #Create output
        if [ \${FailedHEEA} == "Yes" ];then
            > HEEA.BAD.CHECK
            tar cf - HEEA.BAD.CHECK  | \
                pigz -9 -p ${task.cpus} > ./Final/HEEAlPack.tar.gz
        else
            > HEEA.GOOD.CHECK
            tar cf - ./Round*/ HEEA.GOOD.CHECK ./Final/ | \
                pigz -9 -p ${task.cpus} > ./Final/HEEA.tar.gz
        fi

        #Make final check
        > ./Final/GoodCHECK
    fi

    #Link results to oldDir
    cd \${oldWorkDir}
    ln -s \${NewPath}/Final/HEEAPack.tar.gz ./
    > FINAL_HEEA.CHECK
    """
}