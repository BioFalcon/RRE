process RRE_CentralExtension{
    label "RRE_Extension"

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
    file(ExtractCoordScript)

    output:
    tuple env(Repeat),file("Central.FINISH.CHECK"), emit:AllSides

    script:
    """
    #Load ID of the repeat into a variable
    Repeat=\$( cat ${Consensus} )
    FullID=\$( grep "\${Repeat}#" ${ConsensusFasta} | sed 's/>//;s/ .*//' )

    #Save nf directory into a variable
    oldWorkDir=\$(pwd)

    #Move to new WorkDir 
    cd ${outDir}/WorkDir
    mkdir -p \${Repeat}/Central
    cd \${Repeat}/Central

    #Backup new path
    NewPath=\$(pwd)
    #Establish how many sequences will be used
    SampleSeqs=${params.SampleSeqs}
    CoverageKeep=0.25

    #Initialize variables
    ExtensionSize=${params.extension}
    Increase=\$(echo "\${ExtensionSize}" | awk '{\$0=\$0*0.25;printf("%d\\n",\$0+=\$0<0?0:0.9)}' )

    #Family variables
    FamilyNum=00
    NumFamilies=1
    
    #Check if Staging dir exist
    if [ ! -d Staging ] || [ ! -f ./Staging/GoodCHECK ] ;then
        mkdir -p Staging
        ln -s -f \${oldWorkDir}/* ./Staging/

        #Extract the hmmer entries for that repeat
        awk -v Repeat="\${Repeat}" '{OFS="\\t"} \$3==Repeat{print}' ./Staging/${HMMROut} > ./Staging/HMMROut_Subset.txt

        #Make check for Staging
        > ./Staging/GoodCHECK
    fi

    #Do check before 
    if [ -f Central.GOOD.CHECK ];then
        FamilNum=\$( echo NumFamilies | awk '{printf "%02d\\n", \$1+1}' )
    fi
    
    while [[ \${FamilyNum} -lt \${NumFamilies} ]];do
        RoundCentral=01
        GoodCheck="Yes"
        SkipRounds=No

        #Init checks
        FailedCentral=No
        FailedLeft=No
        FailedRight=No
        FinalFiles=No

        #Create Family dir
        mkdir -p ./Family_\${FamilyNum}

        #Make the first round
        if [ ! -d ./Family_\${FamilyNum}/Round_00 ] || [ ! -f ./Family_\${FamilyNum}/Round_00/GoodCHECK ]; then
            #Create dir
            mkdir -p ./Family_\${FamilyNum}/Round_00
            
            #Select which sequences will be used 
            sort -k14,14nr ./Staging/HMMROut_Subset.txt | \\
            cut -f1,9,10,12,14 | \\
            awk '{OFS="\\t"}{ if (\$4 == "+"){print \$1,\$2,\$3,".",\$5,\$4}else{print \$1,\$3,\$2,".",\$5,\$4}}' | \\
            sort -k5,5nr > ./Family_\${FamilyNum}/Round_00/02_CurrentConsensi.Round00.Extended.bed 

            #Extract previous consensus
            samtools faidx ./Staging/${ConsensusFasta} \${FullID} > ./Family_\${FamilyNum}/Round_00/07_CurrentConsensi.Round00.Extended.Curated.Consensus.fa

            #Make check
            > ./Family_\${FamilyNum}/Round_00/GoodCHECK
        fi

        #Check if other Rounds exist and backup if needed
        if [ \$(ls -d ./Family_\${FamilyNum}/Round_*| wc -l) -gt 1 ];then
            #Check if Final dir is there
            if [ -d ./Family_\${FamilyNum}/Final ] ;then
                #Set variable to not do rounds
                SkipRounds=Yes
                #Check if Final is OK
                if [ -f ./Family_\${FamilyNum}/GoodCHECK ];then
                    FinalFiles=Yes
                fi
            else
                #Get which is the last round, and back it up
                RoundCentral=\$(ls -d ./Family_\${FamilyNum}/Round_*| tail -n 1| sed 's/.*Round_//')

                #Check if it is complete and increase round
                if [ -f ./Family_\${FamilyNum}/Round_\${RoundCentral}/GoodCHECK ] && [ \${RoundCentral} -lt ${params.maxCentralRounds} ];then
                    RoundCentral=\$(echo \${RoundCentral} | awk '{printf "%02d\\n", \$1+1}')
                elif [ -f ./Family_\${FamilyNum}/Round_\${RoundCentral}/GoodCHECK ] && [ \${RoundCentral} -eq ${params.maxCentralRounds} ];then
                    SkipRounds=Yes
                    RoundCentral=\$(echo \${RoundCentral} | awk '{printf "%02d\\n", \$1+1}')
                else
                    #Make backup of failed round 
                    BackupNum=01
                    Back=True
                    while [ \${Back} == "True" ];do
                        if [ -d ./Family_\${FamilyNum}/BackUp_Round\${RoundCentral}_\${BackupNum} ];then
                            BackupNum=\$(echo \${BackupNum} | awk '{printf "%02d\\n", \$1+1}')
                        else
                            mv ./Family_\${FamilyNum}/Round_\${RoundCentral} ./Family_\${FamilyNum}/BackUp_Round\${RoundCentral}_\${BackupNum}
                            Back=False
                        fi
                    done
                fi
            fi
        fi

        while [[ \${RoundCentral} -le ${params.maxCentralRounds} ]] && [ \${SkipRounds} == "No" ]; do
            #Make dir for this round
            mkdir -p ./Family_\${FamilyNum}/Round_\${RoundCentral}

            #Set previous round
            PrevRoundCentral=\$(echo \${RoundCentral} | awk '{printf "%02d\\n", \$1-1}')

            #Extract coordinates in sense
            ln -s -f \$(realpath ./Family_\${FamilyNum}/Round_\${PrevRoundCentral}/02_CurrentConsensi.Round\${PrevRoundCentral}.Extended.bed) ./Family_\${FamilyNum}/Round_\${RoundCentral}/00_CurrentConsensi.Round\${RoundCentral}.coord.bed
        
            #Extend coordinates
            bedtools slop \\
                -b \${ExtensionSize} \\
                -i ./Family_\${FamilyNum}/Round_\${RoundCentral}/00_CurrentConsensi.Round\${RoundCentral}.coord.bed \\
                -g <(cut -f1,2 ./Staging/${GenomeIndex}) \\
                > ./Family_\${FamilyNum}/Round_\${RoundCentral}/01_CurrentConsensi.Round\${RoundCentral}.Extended.bed

            #Select the sequences to be used
            CounterSeqs=0
            FileLine=1
            TotalLines=\$(wc -l < ./Family_\${FamilyNum}/Round_\${RoundCentral}/01_CurrentConsensi.Round\${RoundCentral}.Extended.bed)
            > ./Family_\${FamilyNum}/Round_\${RoundCentral}/02_CurrentConsensi.Round\${RoundCentral}.Extended.bed

            while [ \$CounterSeqs -lt \${SampleSeqs} ];do
                #Extract current line
                head -n \${FileLine} ./Family_\${FamilyNum}/Round_\${RoundCentral}/01_CurrentConsensi.Round\${RoundCentral}.Extended.bed | tail -n 1  > ./Family_\${FamilyNum}/Round_\${RoundCentral}/Temp.bed

                #Check if there is no overlap
                if [[ \$(bedtools intersect -a ./Family_\${FamilyNum}/Round_\${RoundCentral}/Temp.bed -b ./Family_\${FamilyNum}/Round_\${RoundCentral}/02_CurrentConsensi.Round\${RoundCentral}.Extended.bed | wc -l) -eq 0 ]] ;then
                    cat ./Family_\${FamilyNum}/Round_\${RoundCentral}/Temp.bed >> ./Family_\${FamilyNum}/Round_\${RoundCentral}/02_CurrentConsensi.Round\${RoundCentral}.Extended.bed
                    CounterSeqs=\$((\${CounterSeqs}+1))
                fi

                #Increase FileLine
                FileLine=\$((\${FileLine}+1))

                if [[ \${FileLine} -gt \${TotalLines} ]];then
                    break
                fi
            done

            #Check if it has enough sequences to continue
            if [ \$(wc -l < ./Family_\${FamilyNum}/Round_\${RoundCentral}/02_CurrentConsensi.Round\${RoundCentral}.Extended.bed ) -lt 10 ] ;then
                break
            fi

            #Get fasta sequences
            bedtools getfasta -fi ./Staging/${Genome} -fo ./Family_\${FamilyNum}/Round_\${RoundCentral}/03_CurrentConsensi.Round\${RoundCentral}.Extended.fa -bed ./Family_\${FamilyNum}/Round_\${RoundCentral}/02_CurrentConsensi.Round\${RoundCentral}.Extended.bed -s

            #Set sequence filter
            ##At least 20 percent of the sequences used
            SequenceFilter=\$(wc -l < ./Family_\${FamilyNum}/Round_\${RoundCentral}/02_CurrentConsensi.Round\${RoundCentral}.Extended.bed | awk '{\$0=\$0*0.2;printf("%d\\n",\$0+=\$0<0?0:0.9)}')

            #Make alignment with MAFFT
            mafft \\
                --genafpair \\
                --maxiterate 1000 \\
                --thread ${task.cpus} \\
                ./Family_\${FamilyNum}/Round_\${RoundCentral}/03_CurrentConsensi.Round\${RoundCentral}.Extended.fa \\
            > ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_CurrentConsensi.Round\${RoundCentral}.Extended.aln

            #Curate alignment
            python3 ./Staging/${CurateScript} \\
                --input ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_CurrentConsensi.Round\${RoundCentral}.Extended.aln \\
                --output ./Family_\${FamilyNum}/Round_\${RoundCentral}/05_CurrentConsensi.Round\${RoundCentral}.Extended.Curated \\
                --StepWindow 2 \\
                --Plots \\
                --MinCoverage 15 \\
                --NoiseThresh 0.05

            cp ./Family_\${FamilyNum}/Round_\${RoundCentral}/05_CurrentConsensi.Round\${RoundCentral}.Extended.Curated__Peak0.aln.fa \\
                ./Family_\${FamilyNum}/Round_\${RoundCentral}/07_CurrentConsensi.Round\${RoundCentral}.Extended.Curated.aln.fa

            #Replace Consensus with hmmer one
            hmmbuild \\
                --dna \\
                --informat afa \\
                --seed 1992 \\
                --cpu ${task.cpus} \\
                --symfrac 0.20 \\
                --fragthresh 1.0 \\
                --wnone \\
                -n \${Repeat} \\
                ./Family_\${FamilyNum}/Round_\${RoundCentral}/07_CurrentConsensi.Round\${RoundCentral}.Extended.Curated.hmm \\
                ./Family_\${FamilyNum}/Round_\${RoundCentral}/07_CurrentConsensi.Round\${RoundCentral}.Extended.Curated.aln.fa

            hmmemit \\
                -c \\
                ./Family_\${FamilyNum}/Round_\${RoundCentral}/07_CurrentConsensi.Round\${RoundCentral}.Extended.Curated.hmm \\
                > ./Family_\${FamilyNum}/Round_\${RoundCentral}/07_CurrentConsensi.Round\${RoundCentral}.Extended.Curated.Consensus.fa
                    
            #Check if it should be terminated
            PrevRoundCentral=\$(echo \${RoundCentral} | awk '{printf "%02d\\n", \$1-1}')
            Length1=\$( grep -v "^>" ./Family_\${FamilyNum}/Round_\${PrevRoundCentral}/07_CurrentConsensi.Round\${PrevRoundCentral}.Extended.Curated.Consensus.fa | sed 's/N//g'| tr '\\n' ' '| sed 's/ \$/\\n/;s/ //g' | wc -m )
            Length2=\$( echo "\$(grep -v '^>' ./Family_\${FamilyNum}/Round_\${RoundCentral}/07_CurrentConsensi.Round\${RoundCentral}.Extended.Curated.Consensus.fa | sed 's/N//g'| tr '\\n' ' '| sed 's/ \$/\\n/;s/ //g' | wc -m ) - \${Increase}"| bc)

            #If it made it this far, make check 
            > ./Family_\${FamilyNum}/Round_\${RoundCentral}/GoodCHECK

            if (( \$(echo " \$Length2 < \$Length1" | bc -l) ));then 
                break
            else
                #Make check for this round
                #> ./Round_\${RoundCentral}/GoodCHECK
                RoundCentral=\$(echo \${RoundCentral} | awk '{printf "%02d\\n", \$1+1}')
            fi
        done

        #############################
        if [ \${FinalFiles} == "No" ];then
            #Make Final directory
            mkdir -p ./Family_\${FamilyNum}/Final

            #Determine last OK round
            LastSuccRound=\$(ls ./Family_\${FamilyNum}/Round_*/GoodCHECK| tail -1| sed 's/\\/GoodCHECK//;s/.*Round_//')
            PrevSucRound=\$( echo \${LastSuccRound} | awk '{printf "%02d\\n", \$1-1}' )

            if [ \${LastSuccRound} == '00' ];then
                echo ">\${Repeat}#___Failed" > ./Family_\${FamilyNum}/FinalConsensi.aln.fa
                echo "A" >> ./Family_\${FamilyNum}/Final/FinalConsensi.aln.fa
                echo ">\${Repeat}#___Failed" > ./Family_\${FamilyNum}/Final/FinalConsensi.consensus.fa
                echo "A" >> ./Family_\${FamilyNum}/Final/FinalConsensi.consensus.fa
                FailedCentral=Yes
            else
                cp ./Family_\${FamilyNum}/Round_\${LastSuccRound}/07_CurrentConsensi.Round\${LastSuccRound}.Extended.Curated.aln.fa \\
                    ./Family_\${FamilyNum}/Final/FinalConsensi.aln.fa
                cp ./Family_\${FamilyNum}/Round_\${LastSuccRound}/07_CurrentConsensi.Round\${LastSuccRound}.Extended.Curated.Consensus.fa \\
                    ./Family_\${FamilyNum}/Final/FinalConsensi.consensus.fa
            fi

            #Copy coords that were used
            cp ./Family_\${FamilyNum}/Round_\${LastSuccRound}/02_CurrentConsensi.Round\${LastSuccRound}.Extended.bed \\
                ./Family_\${FamilyNum}/Final/FinalConsensi.Extension.Coord.bed

            #Extract each side of the consensus and decide if should still extend
            if [ \${FailedCentral} == "No" ];then
                #Make MAFFT alignment
                mafft \\
                    --localpair \\
                    --maxiterate 1000 \\
                    --thread ${task.cpus} \\
                    <( cat ./Family_\${FamilyNum}/Round_\${PrevSucRound}/07_CurrentConsensi.Round\${PrevSucRound}.Extended.Curated.Consensus.fa \\
                        ./Family_\${FamilyNum}/Round_\${LastSuccRound}/07_CurrentConsensi.Round\${LastSuccRound}.Extended.Curated.Consensus.fa ) \\
                    > ./Family_\${FamilyNum}/Final/FinalSideCheck.aln.fa

                #Extract coordinates from alignment
                ./Staging/${ExtractCoordScript} \\
                    --input ./Family_\${FamilyNum}/Final/FinalSideCheck.aln.fa \\
                    --output ./Family_\${FamilyNum}/Final/FinalConsensi.SidesCheck.Coords \\
                    --Side Central

                ###Get Length of extension into variables
                LengthExtension=\$(sed '/>/d' ./Family_\${FamilyNum}/Final/FinalConsensi.consensus.fa| tr '\\n' ' '| sed 's/ //g'| wc -m)
                LeftSide_Size=\$(head -n 1 ./Family_\${FamilyNum}/Final/FinalConsensi.SidesCheck.Coords | awk '{print \$2 - \$1}')
                RightSide_Size=\$(tail -n 1 ./Family_\${FamilyNum}/Final/FinalConsensi.SidesCheck.Coords | awk '{print \$2 - \$1}')

                ##CheckSides
                ###Left
                if [ \${LeftSide_Size} -lt \${Increase} ]; then
                    FailedLeft=Yes
                else
                    head -n 1 ./Family_\${FamilyNum}/Final/FinalConsensi.SidesCheck.Coords > ./Family_\${FamilyNum}/Final/Final.Consensus.LeftSide.Coord
                fi

                ###Right
                if [ \${RightSide_Size} -lt \${Increase} ]; then
                    FailedRight=Yes
                else
                    tail -n 1 ./Family_\${FamilyNum}/Final/FinalConsensi.SidesCheck.Coords > ./Family_\${FamilyNum}/Final/Final.Consensus.RightSide.Coord
                fi

            fi

            #Pack-up results
            mkdir -p ./Family_\${FamilyNum}/Final/LeftPackCentral
            mkdir -p ./Family_\${FamilyNum}/Final/RightPackCentral
            if [ \${FailedCentral} == "Yes" ];then
                > ./Family_\${FamilyNum}/Final/LeftPackCentral/LEFT.BAD.CHECK
                > ./Family_\${FamilyNum}/Final/RightPackCentral/RIGHT.BAD.CHECK
            else
                if [ \${FailedLeft} == "Yes" ];then
                    > ./Family_\${FamilyNum}/Final/LeftPackCentral/LEFT.BAD.CHECK
                else 
                    cp ./Family_\${FamilyNum}/Final/FinalConsensi.aln.fa \\
                        ./Family_\${FamilyNum}/Final/LeftPackCentral/
                    cp ./Family_\${FamilyNum}/Final/FinalConsensi.consensus.fa \\
                        ./Family_\${FamilyNum}/Final/LeftPackCentral/
                    cp ./Family_\${FamilyNum}/Final/Final.Consensus.LeftSide.Coord \\
                        ./Family_\${FamilyNum}/Final/LeftPackCentral/
                    > ./Family_\${FamilyNum}/Final/LeftPackCentral/LEFT.GOOD.CHECK
                fi
                if [ \${FailedRight} == "Yes" ];then
                    > ./Family_\${FamilyNum}/Final/RightPackCentral/RIGHT.BAD.CHECK
                else
                    cp ./Family_\${FamilyNum}/Final/FinalConsensi.aln.fa \\
                        ./Family_\${FamilyNum}/Final/RightPackCentral/
                    cp ./Family_\${FamilyNum}/Final/FinalConsensi.consensus.fa \\
                        ./Family_\${FamilyNum}/Final/RightPackCentral/
                    cp ./Family_\${FamilyNum}/Final/Final.Consensus.RightSide.Coord \\
                        ./Family_\${FamilyNum}/Final/RightPackCentral/
                    > ./Family_\${FamilyNum}/Final/RightPackCentral/RIGHT.GOOD.CHECK
                fi
            fi

            #Compress files
            tar cf - ./Family_\${FamilyNum}/Final/LeftPackCentral  | \\
                pigz -9 -p ${task.cpus} > ./Family_\${FamilyNum}/Final/LeftPack.Central.tar.gz
            tar cf - ./Family_\${FamilyNum}/Final/RightPackCentral | \\
                pigz -9 -p ${task.cpus} > ./Family_\${FamilyNum}/Final/RightPack.Central.tar.gz

            #Create output for central piece
            if [ \${FailedCentral} == "Yes" ];then
                > ./Family_\${FamilyNum}/Central.BAD.CHECK
                tar cf - Central.BAD.CHECK  | \\
                    pigz -9 -p ${task.cpus} > ./Family_\${FamilyNum}/Final/CentralPack.tar.gz
            else
                > ./Family_\${FamilyNum}/Central.GOOD.CHECK
                tar cf - ./Family_\${FamilyNum}/Round*/ Central.GOOD.CHECK ./Family_\${FamilyNum}/Final/ | \\
                    pigz -9 -p ${task.cpus} > ./Family_\${FamilyNum}/Final/CentralPack.tar.gz
            fi

            #Make final check
            > ./Family_\${FamilyNum}/Final/GoodCHECK
        fi

        #Update variables
        FamilyNum=\$( echo \${FamilyNum} | awk '{printf "%02d\\n", \$1+1}' )

        NumFamilies=\$(ls -d ./Family_* | wc -l)
    done

    #Make Check for finished step
    > ./Central.FINISH.CHECK

    #Link results to oldDir
    ln -s \$( realpath ./Central.FINISH.CHECK ) \${oldWorkDir}
    cd \${oldWorkDir}
    """
}