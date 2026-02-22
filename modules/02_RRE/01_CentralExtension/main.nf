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

            #Make alignment
            mafft \\
                --localpair \\
                --maxiterate 1000 \\
                --thread ${task.cpus} \\
                ./Family_\${FamilyNum}/Round_\${RoundCentral}/03_CurrentConsensi.Round\${RoundCentral}.Extended.fa \\
                > ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_CurrentConsensi.Round\${RoundCentral}.Extended.aln \\
                2> ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_CurrentConsensi.Round\${RoundCentral}.Extended.MAFFTlog

            #Extract and uppercase 
            seqkit seq \\
                -u \\
                ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_CurrentConsensi.Round\${RoundCentral}.Extended.aln \\
                > 04_TEMP.aln
            mv 04_TEMP.aln ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_CurrentConsensi.Round\${RoundCentral}.Extended.aln

            #Check if there are other families
            ./Staging/${OutlierScript} \\
                --input ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_CurrentConsensi.Round\${RoundCentral}.Extended.aln \\
                --output ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_CurrentConsensi.Round\${RoundCentral}.Extended.Fams

            #Check if there are multiple families
            if [  \$(ls ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_CurrentConsensi.Round\${RoundCentral}.Extended.Fams.Cluster*fa.aln | wc -l ) -gt 1 ];then
                #Only take those with more than N sequences
                CurrFamily=00
                for FamilyAln in \$(ls ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_CurrentConsensi.Round\${RoundCentral}.Extended.Fams.Cluster*fa.aln );do
                        if [[ \$(grep -c ">" \${FamilyAln}) -ge 7 ]];then
                            cp \${FamilyAln} ./Family_\${FamilyNum}/Round_\${RoundCentral}/4_00_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}.aln.fa
                        fi
                        CurrFamily=\$( echo \${CurrFamily} | awk '{printf "%02d\\n", \$1+1}' )
                done

                #Break if no families have more than n sequences
                if [ \$( grep -c ">" ./Family_\${FamilyNum}/Round_\${RoundCentral}/4_00_CurrentConsensi.Round\${RoundCentral}.Family*.aln.fa | sed 's/.*://' | awk '\$1>10{print}' | wc -l ) -lt 1 ];then
                    FailedCentral=Yes
                    break
                fi

                #Loop through families
                for FamilyAln in \$( ls ./Family_\${FamilyNum}/Round_\${RoundCentral}/4_00_CurrentConsensi.Round\${RoundCentral}.Family*.aln.fa );do
                    #Get family number
                    CurrFamily=\$( echo \${FamilyAln} | \\
                            sed 's/.*Family//;s/.aln.fa//' )

                    #Clean alignment
                    ./Staging/${CurateScript} \\
                    --input \${FamilyAln} \\
                    --output ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_01_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily} \\
                    --SeqID \${Repeat}#___Round\${RoundCentral} \\
                    --MinCoverage 1 \\
                    --Plots

                    #Select best peak
                    PeakNumLoop=0

                    for PeakAln in \$( ls ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_01_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}__Peak*aln.fa );do
                        if [[ \$(seqkit stats -T \${PeakAln} | cut -f8 | tail -1) -gt 0 ]];then 
                            #Create hmmer model
                            hmmbuild \\
                            --symfrac 0 \\
                            --dna \\
                            --informat afa \\
                            --seed 1992 \\
                            --fragthresh 1.0 \\
                            --wnone \\
                            -n Peak\${PeakNumLoop} \\
                            ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_02_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}_Peak\${PeakNumLoop}.hmm \\
                            \${PeakAln}

                            #Search for the model in the genome
                            nhmmer \\
                            --cpu ${task.cpus} \\
                            --tblout ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_03_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}_Peak\${PeakNumLoop}.tblout \\
                            --notextw \\
                            --noali \\
                            --seed 1992 \\
                            -E 1e-5 \\
                            ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_02_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}_Peak\${PeakNumLoop}.hmm \\
                            ./Staging/HMMDB_Genome > /dev/null
                        
                            #Cleanup and remove high bias
                            sed '/#/d; s/\\s\\{1,\\}/\\t/g' ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_03_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}_Peak\${PeakNumLoop}.tblout | \\
                            awk '{OFS="\\t"} (\$15/\$14) <= 0.20{print}' > ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_03_Temp.tblout 
                            mv ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_03_Temp.tblout \\
                                ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_03_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}_Peak\${PeakNumLoop}.tblout
                        else
                            > ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_03_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}_Peak\${PeakNumLoop}.tblout
                            > ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_02_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}_Peak\${PeakNumLoop}.hmm
                        fi

                        #Increase peak number
                        PeakNumLoop=\$((\${PeakNumLoop}+1))
                    done

                    #Select Peak with best hits
                    PeakNumLoop=0
                    HighestMedian=0
                    SelectedPeak=0
                    for PeakHMMOut in \$(ls ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_03_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}_Peak*.tblout );do
                        #Calculate median of bitscores
                        MedianScores=\$( sed 's/\\s\\{1,\\}/\\t/g' \${PeakHMMOut} | \\
                                        cut -f14 | \\
                                        sort -n | \\
                                        awk '{a[i++]=\$1} END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]}' )

                            #Check if median is higher than the highest median
                        if (( \$(echo "\${MedianScores} > \${HighestMedian}" | bc -l) )) ;then
                            HighestMedian=\${MedianScores}
                            SelectedPeak=\${PeakNumLoop}
                        fi
                        PeakNumLoop=\$((\${PeakNumLoop}+1))
                    done

                    #Select the best peak and hits
                    cp ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_01_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}__Peak\${SelectedPeak}.aln.fa \\
                    ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_04_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}.aln.fa

                    cp ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_03_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}_Peak\${SelectedPeak}.tblout \\
                    ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_04_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}.tblout

                    cp ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_01_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}.Consensus__Peak\${SelectedPeak}.fa \\
                    ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_04_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}.Consensus.fa

                    #Select sequences that cover N% of the model
                    ModelLen=\$(seqkit stats -T ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_04_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}.aln.fa | \\
                            cut -f7 | \\
                            tail -1 | \\
                            sed 's/\\..*//' )

                    bedtools intersect \\
                    -a <(echo -e "Target\\t0\\t\${ModelLen}") \\
                    -b <( sed '/#/d' ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_04_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}.tblout |  awk '{OFS="\\t"}{print "Target",\$5,\$6,\$1,\$9,\$10,\$12,\$14}' ) \\
                    -wb \\
                    -f \${CoverageKeep} | \\
                    awk '{OFS="\\t"}{ 
                    if( \$10 == "+"){print \$7,\$8,\$9 ,".",\$11,\$10}else{print \$7,\$9,\$8,".",\$11,\$10}
                    }' > ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_05_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}.bed

                    #Select sequences to use
                    CounterSeqs=0
                    FileLine=1
                    TotalLines=\$(wc -l < ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_05_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}.bed )
                    > ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_06_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}.bed

                    while [[ \${CounterSeqs} -lt \${SampleSeqs} ]];do
                        #Extract current line
                        head -n \${FileLine} ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_05_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}.bed | tail -n 1  > ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_Temp.bed

                            #Check if there is no overlap
                        if [[ \$(bedtools intersect -a ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_Temp.bed -b ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_06_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}.bed | wc -l ) -eq 0 ]];then
                            cat ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_Temp.bed >>  ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_06_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}.bed
                            CounterSeqs=\$((\${CounterSeqs}+1))
                        fi

                        #Increase FileLine
                        FileLine=\$((\${FileLine}+1))

                        if [[ \${FileLine} -gt \${TotalLines} ]];then
                            break
                        fi
                    done

                    #Extract sequences
                    bedtools getfasta \\
                    -fi  ./Staging/${Genome} \\
                    -fo  ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_07_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}.fa \\
                    -bed ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_06_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}.bed \\
                    -s

                    #Add ID about Round
                    awk -v Round=\${RoundCentral} \\
                    '{ 
                            if(\$1 ~ "^>"){
                            gsub("\$","___Central"Round,\$1);
                            print
                            } else{print}
                    }' \\
                    ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_07_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}.fa \\
                    > ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_Temp.fa
                    
                    mv ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_Temp.fa \\
                        ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_07_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}.fa

                    #Make alignment for this round
                    mafft \\
                    --localpair \\
                    --maxiterate 1000 \\
                    --thread ${task.cpus} \\
                    ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_07_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}.fa \\
                    > ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_08_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}.aln

                    #Make all sequences uppercase
                    seqkit seq \\
                    -u \\
                    ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_08_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}.aln \\
                    > ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_08_Temp.aln

                    mv ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_08_Temp.aln \\
                    ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_08_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}.aln

                    #Curate alignment and make consensus
                    ./Staging/${CurateScript} \\
                    --input ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_08_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}.aln \\
                    --output ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_09_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily} \\
                    --SeqID \${Repeat}#___RoundCentral\${RoundCentral} \\
                    --NoiseThresh ${params.noiseThreshold} \\
                    --Plots

                    #Search each peak in genome
                    PeakNumLoop=0

                    for PeakAln in \$(ls ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_09_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}__Peak*aln.fa );do
                        if [[ \$(seqkit stats -T \${PeakAln} | cut -f8 | tail -1) -gt 0 ]];then 
                            #Create hmmer model
                            hmmbuild \\
                            --symfrac 0 \\
                            --dna \\
                            --informat afa \\
                            --seed 1992 \\
                            --fragthresh 1.0 \\
                            --wnone \\
                            -n Peak\${PeakNumLoop} \\
                            ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_10_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}_Peak\${PeakNumLoop}.hmm \\
                            \${PeakAln}

                            #Search for the model in the genome
                            nhmmer \\
                            --cpu ${task.cpus} \\
                            --tblout ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_11_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}_Peak\${PeakNumLoop}.tblout \\
                            --notextw \\
                            --noali \\
                            --seed 1992 \\
                            -E 1e-5 \\
                            ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_10_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}_Peak\${PeakNumLoop}.hmm \\
                            ./Staging/HMMDB_Genome > /dev/null

                            #Cleanup and remove high bias
                            sed '/#/d; s/\\s\\{1,\\}/\\t/g' ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_11_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}_Peak\${PeakNumLoop}.tblout | \\
                            awk '{OFS="\\t"} (\$15/\$14) <= 0.20{print}' > ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_11_Temp.tblout 
                            mv ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_11_Temp.tblout ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_11_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}_Peak\${PeakNumLoop}.tblout
                        else
                            > ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_11_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}_Peak\${PeakNumLoop}.tblout
                            > ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_10_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}_Peak\${PeakNumLoop}.hmm 
                        fi

                        #Select Peak with best hits
                        PeakNumLoop=0
                        HighestMedian=0
                        SelectedPeak=0

                        for PeakHMMOut in \$(ls ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_11_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}_Peak*.tblout );do
                            #Calculate median of bitscores
                            MedianScores=\$( sed 's/\\s\\{1,\\}/\\t/g' \${PeakHMMOut} | \\
                                            cut -f14 | \\
                                            sort -n | \\
                                            awk '{a[i++]=\$1} END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]}' )

                            #Check if median is higher than the highest median
                            if (( \$(echo "\${MedianScores} > \${HighestMedian}" | bc -l) )) ;then
                                HighestMedian=\${MedianScores}
                                SelectedPeak=\${PeakNumLoop}
                            fi
                            PeakNumLoop=\$((\${PeakNumLoop}+1))
                        done
                        
                        #Select the best peak and hits
                        cp ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_09_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}__Peak\${SelectedPeak}.aln.fa \\
                        ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_12_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}.aln.fa

                        cp ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_11_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}_Peak\${SelectedPeak}.tblout \\
                        ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_12_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}.tblout

                        cp ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_09_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}.Consensus__Peak\${SelectedPeak}.fa \\
                        ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_12_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}.Consensus.fa

                    done

                    #Select those families that have homology with the previous round
                    ##Only do if round is greater than 1
                    if [[ \${RoundCentral} -gt 1 ]];then
                        PrevRound=\$(echo \${RoundCentral} | awk '{printf "%02d\\n", \$1-1}')

                        #Make a db from previous consensus
                        makeblastdb \\
                                -in ./Family_\${FamilyNum}/Round_\${RoundCentral}/07_CurrentConsensi.Round\${RoundCentral}.Extended.Curated.Consensus.fa \\
                                -dbtype nucl \\
                                -out ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_13_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}.db
                        
                        #Only do blast if there is sequence left
                        if [[ \$( seqkit stat -T ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_12_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}.Consensus.fa | cut -f7| tail -1| sed 's/\\..*//' ) -gt 0 ]];then
                            #Perform blast
                            blastn \\
                            -query ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_12_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}.Consensus.fa \\
                            -db ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_13_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}.db \\
                            -outfmt 6 \\
                            -out ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_13_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}.blast \\
                            -num_threads ${task.cpus}

                            #if there are hits, then keep the family
                            if [[ \$(wc -l < ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_13_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}.blast) -gt 0 ]];then
                                cp ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_12_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}.aln.fa ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_14_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}.aln.fa
                            fi
                        fi
                    else
                        cp ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_12_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}.aln.fa ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_14_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}.aln.fa
                    fi

                    #Increase family number
                    CurrFamily_Num=\$( echo \${CurrFamily} | awk '{printf "%02d\\n", \$1+1}' )
                done

                if [ \$(ls ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_14_CurrentConsensi.Round\${RoundCentral}.Family*.aln.fa | wc -l ) -gt 0 ];then
                    #Loop through families
                    for CurrFamily in \$(ls ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_14_CurrentConsensi.Round\${RoundCentral}.Family*.aln.fa | sed 's/.*Family//;s/.aln.fa//');do
                        #Do all the steps for this one before sorting
                        ##Cleanup alignment
                        ./Staging/${CurateScript} \\
                            --input ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_14_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}.aln.fa \\
                            --output ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_15_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily} \\
                            --SeqID \${Repeat}#___RoundC\${RoundCentral} \\
                            --NoiseThresh ${params.noiseThreshold} \\
                            --Plots
                        
                        ##Select the best peak
                        PeakNumLoop=0
                        if [[ \$(ls ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_15_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}__Peak*aln.fa | wc -l ) -gt 1 ]];then
                            for PeakAln in \$(ls ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_15_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}__Peak*aln.fa );do
                                if [[ \$(seqkit stats -T \${PeakAln} | cut -f8 | tail -1) -gt 0 ]];then 
                                    #Create hmmer model
                                    hmmbuild \\
                                    --symfrac 0 \\
                                    --dna \\
                                    --informat afa \\
                                    --seed 1992 \\
                                    --fragthresh 1.0 \\
                                    --wnone \\
                                    -n Peak\${PeakNumLoop} \\
                                    ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_15_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}_Peak\${PeakNumLoop}.hmm \\
                                    \${PeakAln}

                                    #Search for the model in the genome
                                    nhmmer \\
                                    --cpu ${task.cpus} \\
                                    --tblout ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_15_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}_Peak\${PeakNumLoop}.tblout \\
                                    --notextw \\
                                    --noali \\
                                    --seed 1992 \\
                                    -E 1e-5 \\
                                    ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_15_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}_Peak\${PeakNumLoop}.hmm \\
                                    ./Staging/HMMDB_Genome > /dev/null

                                    #Cleanup and remove high bias
                                    sed '/#/d; s/\\s\\{1,\\}/\\t/g' ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_15_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}_Peak\${PeakNumLoop}.tblout | \\
                                    awk '{OFS="\\t"} (\$15/\$14) <= 0.20{print}' \\
                                    > ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_15_Temp.tblout 
                                    
                                    mv ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_15_Temp.tblout \\
                                    ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_15_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}_Peak\${PeakNumLoop}.tblout
                                else
                                    > ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_15_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}_Peak\${PeakNumLoop}.tblout
                                    > ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_15_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}_Peak\${PeakNumLoop}.hmm
                                fi
                                #Increase peak number
                                PeakNumLoop=\$((\${PeakNumLoop}+1))
                            done

                            #Select Peak with best hits
                            PeakNumLoop=0
                            HighestMedian=0
                            SelectedPeak=0

                            for PeakHMMOut in \$(ls ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_15_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}_Peak*.tblout );do
                                #Calculate median of bitscores
                                MedianScores=\$( sed 's/\\s\\{1,\\}/\\t/g' \${PeakHMMOut} | \\
                                                cut -f14 | \\
                                                sort -n | \\
                                                awk '{a[i++]=\$1} END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]}' )

                                #Check if median is higher than the highest median
                                if (( \$(echo "\${MedianScores} > \${HighestMedian}" | bc -l) )) ;then
                                    HighestMedian=\${MedianScores}
                                    SelectedPeak=\${PeakNumLoop}
                                fi
                                PeakNumLoop=\$((\${PeakNumLoop}+1))
                            done
                            
                            #Select the best peak and hits
                            cp ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_15_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}__Peak\${SelectedPeak}.aln.fa \\
                            ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_16_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}.aln.fa

                            cp ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_15_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}_Peak\${SelectedPeak}.tblout \\
                            ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_16_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}.tblout

                            cp ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_15_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}.Consensus__Peak\${SelectedPeak}.fa \\
                            ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_16_CurrentConsensi..Round\${RoundCentral}.Family\${CurrFamily}.Consensus.fa
                        else
                            cp ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_15_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}__Peak0.aln.fa \\
                            ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_16_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}.aln.fa

                            cp ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_15_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}.Consensus__Peak0.fa \\
                            ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_16_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}.Consensus.fa
                        fi

                        #Only realign if it is not round 00
                        if [[ \${RoundCentral} -eq 01 ]];then
                            cp ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_16_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}.aln.fa \\
                                ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_17_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}.aln
                        else
                            #Get previous round
                            PrevRound=\$(echo \${RoundCentral} | awk '{printf "%02d\\n", \$1-1}')

                            ##Add alignment to previous alignment
                            mafft \\
                                --localpair \\
                                --maxiterate 1000 \\
                                --thread ${task.cpus} \\
                                --lexp -1.5 \\
                                --lop 0.5 \\
                                --add ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_16_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}.aln.fa \\
                                ./Family_\${FamilyNum}/Round_\${PrevRound}/07_CurrentConsensi.Round\${PrevRound}.Extended.Curated.aln.fa \\
                                > ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_17_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}.aln

                        fi

                        #Clean alignment
                        #Make zero consensus and cleaning
                        ./Staging/${CurateScript} \\
                        --input ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_17_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}.aln \\
                        --output ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_18_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}.Curated \\
                        --SeqID \${Repeat}#___RoundC\${RoundCentral} \\
                        --Plots \\
                        --zeroOnly
                        
                        #Remove peak suffix
                        mv ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_18_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}.Curated.Consensus__Peak0.fa \\
                            ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_18_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}.Curated.Consensus.fa
                        mv ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_18_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}.Curated__Peak0.aln.fa \\
                            ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_18_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}.Curated.aln.fa

                        #Check if consensus has grown too much
                        CurrConsensLength=\$( seqkit fx2tab --length --name ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_18_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}.Curated.Consensus.fa | cut -f2 )

                        ##Previous consensus length
                        PrevRound=\$(echo \${RoundCentral} | awk '{printf "%02d\\n", \$1-1}')
                        PrevConsensLength=\$( seqkit fx2tab --length --name ./Family_\${FamilyNum}/Round_\${PrevRound}/07_CurrentConsensi.Round\${PrevRound}.Extended.Curated.Consensus.fa | cut -f2 )

                        ##Make mafft alignment of the previous consensus
                        mafft \\
                            --localpair \\
                            --maxiterate 1000 \\
                            --thread ${task.cpus} \\
                            <( cat ./Family_\${FamilyNum}/Round_\${PrevRound}/07_CurrentConsensi.Round\${PrevRound}.Extended.Curated.Consensus.fa \\
                               ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_18_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}.Curated.Consensus.fa ) \\
                            > ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_19_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}.Curated.aln
                                
                        ./Staging/${ExtractCoordScript} \\
                            --input ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_19_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}.Curated.aln \\
                            --output ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_19_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}.Curated.Coord \\
                            --Side Central

                        
                        > ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_20.GOODCHECK_Family\${CurrFamily}

                        #Make model for next round
                        hmmbuild \\
                        --symfrac 0 \\
                        --dna \\
                        --cpu ${task.cpus} \\
                        --informat afa \\
                        --seed 1992 \\
                        --fragthresh 1.0 \\
                        --wnone \\
                        -n \${Repeat} \\
                        ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_21_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}.Curated.hmm \\
                        ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_18_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFamily}.Curated.aln.fa

                        #Lets sort families into new dirs if they exist
                        FamNum=00
                        if [ \$( ls ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_20.GOODCHECK_Family* | wc -l ) -gt 0 ];then
                                for CurrFam in \$( ls ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_20.GOODCHECK_Family* | sed 's/.*Family//' );do
                                    #First family is reference
                                    if [[ \${FamNum} -eq 0 ]];then
                                        #Fragment alignment 05
                                        cp ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_12_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFam}.aln.fa \\
                                        ./Family_\${FamilyNum}/Round_\${RoundCentral}/05_CurrentConsensi.Round\${RoundCentral}.Extended.Curated.aln.fa
                                        cp ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_12_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFam}.Consensus.fa \\
                                        ./Family_\${FamilyNum}/Round_\${RoundCentral}/05_CurrentConsensi.Round\${RoundCentral}.Extended.Curated.Consensus.fa
                                        
                                        #Add Alignment 06
                                        cp ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_16_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFam}.aln.fa \\
                                        ./Family_\${FamilyNum}/Round_\${RoundCentral}/06_CurrentConsensi.Round\${RoundCentral}.Extended.Curated.aln.fa

                                        #Clean Alignment 07
                                        cp ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_18_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFam}.Curated.aln.fa \\
                                        ./Family_\${FamilyNum}/Round_\${RoundCentral}/07_CurrentConsensi.Round\${RoundCentral}.Extended.Curated.aln.fa
                                        cp ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_18_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFam}.Curated.Consensus.fa \\
                                        ./Family_\${FamilyNum}/Round_\${RoundCentral}/07_CurrentConsensi.Round\${RoundCentral}.Extended.Curated.Consensus.fa

                                        #Coords and model 08
                                        cp ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_19_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFam}.Curated.Coord \\
                                        ./Family_\${FamilyNum}/Round_\${RoundCentral}/08_CurrentConsensi.Round\${RoundCentral}.ExtendedSide.Coord
                                        cp ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_21_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFam}.Curated.hmm \\
                                        ./Family_\${FamilyNum}/Round_\${RoundCentral}/08_CurrentConsensi.Round\${RoundCentral}.Extended.Curated.hmm

                                        #CHECK 
                                        cp ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_20.GOODCHECK_Family\${CurrFam} \\
                                        ./Family_\${FamilyNum}/Round_\${RoundCentral}/GoodCHECK
                                    else
                                        #Make new family directory
                                        ##Last family number
                                        HighestCurFam=\$(ls ./Family_* | sed 's/.*Family_//' | sort -n | tail -1)

                                        ##Increase by one the family number
                                        HighestCurFam=\$( echo \${HighestCurFam} | awk '{printf "%02d\\n", \$1+1}' )

                                        #Make dir
                                        mkdir ./Family_\${HighestCurFam}

                                        #Copy all previous rounds to new family
                                        cp -r ./Family_\${FamilyNum}/Round_* ./Family_\${HighestCurFam}/

                                        #Copy that family into round
                                        ##Fragment alignment 05
                                        cp ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_12_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFam}.aln.fa \\
                                        ./Family_\${HighestCurFam}/Round_\${RoundCentral}/05_CurrentConsensi.Round\${RoundCentral}.Extended.Curated.aln.fa
                                        cp ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_12_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFam}.Consensus.fa \\
                                        ./Family_\${HighestCurFam}/Round_\${RoundCentral}/05_CurrentConsensi.Round\${RoundCentral}.Extended.Curated.Consensus.fa
                                        
                                        ##Add Alignment 06
                                        cp ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_16_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFam}.aln.fa \\
                                        ./Family_\${HighestCurFam}/Round_\${RoundCentral}/06_CurrentConsensi.Round\${RoundCentral}.Extended.Curated.aln.fa

                                        ##Clean Alignment 07
                                        cp ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_18_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFam}.Curated.aln.fa \\
                                        ./Family_\${HighestCurFam}/Round_\${RoundCentral}/07_CurrentConsensi.Round\${RoundCentral}.Extended.Curated.aln.fa
                                        cp ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_18_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFam}.Curated.Consensus.fa \\
                                        ./Family_\${HighestCurFam}/Round_\${RoundCentral}/07_CurrentConsensi.Round\${RoundCentral}.Extended.Curated.Consensus.fa

                                        ##Coords and model 08
                                        cp ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_19_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFam}.Curated.Coord \\
                                        ./Family_\${HighestCurFam}/Round_\${RoundCentral}/08_CurrentConsensi.Round\${RoundCentral}.Extended.Curated.Coords
                                        cp ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_21_CurrentConsensi.Round\${RoundCentral}.Family\${CurrFam}.Curated.hmm \\
                                        ./Family_\${HighestCurFam}/Round_\${RoundCentral}/08_CurrentConsensi.Round\${RoundCentral}.Extended.Curated.hmm

                                        #CHECK 
                                        cp ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_20.GOODCHECK_Family\${CurrFam} \\
                                        ./Family_\${HighestCurFam}/Round_\${RoundCentral}/GoodCHECK
                                    fi
                                    FamNum=\$( echo \${FamNum} | awk '{printf "%02d\\n", \$1+1}' )
                                done
                        else 
                                #Round over
                                BreakNow=True
                        fi
                    done
                else
                    break
                fi  

                #Increase round number
                RoundCentral=\$( echo \${RoundCentral} | awk '{printf "%02d\\n", \$1+1}' )
            else
                #Curate alignment and make consensus
                ./Staging/${CurateScript} \\
                    --input ./Family_\${FamilyNum}/Round_\${RoundCentral}/04_CurrentConsensi.Round\${RoundCentral}.Extended.aln \\
                    --output ./Family_\${FamilyNum}/Round_\${RoundCentral}/05_CurrentConsensi.Round\${RoundCentral}.Extended.Curated \\
                    --SeqID \${Repeat}#___Round\${RoundCentral} \\
                    --NoiseThresh ${params.noiseThreshold} \\
                    --Plots 

                #Select best peak
                PeakNumLoop=0
                for PeakAln in \$(ls ./Family_\${FamilyNum}/Round_\${RoundCentral}/05_CurrentConsensi.Round\${RoundCentral}.Extended.Curated__Peak*aln.fa);do
                    if [[ \$(seqkit stats -T \${PeakAln} | cut -f8 | tail -1) -gt 0 ]];then 
                        #Create hmmer model
                        hmmbuild \\
                        --dna \\
                        --informat afa \\
                        --seed 1992 \\
                        --cpu ${task.cpus} \\
                        --symfrac 0.0 \\
                        --fragthresh 1.0 \\
                        --wnone \\
                        -n Peak\${PeakNumLoop} \\
                        ./Family_\${FamilyNum}/Round_\${RoundCentral}/06_CurrentConsensi__Peak\${PeakNumLoop}.hmm \\
                        \${PeakAln}

                        #Search for the model in the genome
                        nhmmer \\
                        --cpu ${task.cpus} \\
                        --tblout ./Family_\${FamilyNum}/Round_\${RoundCentral}/06_CurrentConsensi__Peak\${PeakNumLoop}.tblout \\
                        --notextw \\
                        --noali \\
                        --seed 1992 \\
                        -E 1e-5 \\
                        ./Family_\${FamilyNum}/Round_\${RoundCentral}/06_CurrentConsensi__Peak\${PeakNumLoop}.hmm \\
                        ./Staging/HMMDB_Genome
                    else
                        > ./Family_\${FamilyNum}/Round_\${RoundCentral}/06_CurrentConsensi__Peak\${PeakNumLoop}.tblout
                        > ./Family_\${FamilyNum}/Round_\${RoundCentral}/06_CurrentConsensi__Peak\${PeakNumLoop}.hmm
                    fi
                        #Increase Peak Number
                        PeakNumLoop=\$((\${PeakNumLoop}+1))
                done

                #Select Peak with best hits
                PeakNumLoop=0
                HighestMedian=0
                SelectedPeak=0
                for PeakHMMOut in \$(ls ./Family_\${FamilyNum}/Round_\${RoundCentral}/06_CurrentConsensi__Peak*.tblout );do
                        #Calculate median of bitscores
                        MedianScores=\$( sed 's/\\s\\{1,\\}/\\t/g' \${PeakHMMOut} | \\
                                        cut -f14 | \\
                                        sort -n | \\
                                        awk '{a[i++]=\$1} END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]}' )

                        #Check if median is higher than the highest median
                        if (( \$(echo "\${MedianScores} > \${HighestMedian}" | bc -l) )) ;then
                                HighestMedian=\${MedianScores}
                                SelectedPeak=\${PeakNumLoop}
                        fi
                        PeakNumLoop=\$((\${PeakNumLoop}+1))
                done

                echo "Selected Peak: \${SelectedPeak}"
                echo "Highest Median: \${HighestMedian}"

                #Copy alignment with the best score
                cp ./Family_\${FamilyNum}/Round_\${RoundCentral}/05_CurrentConsensi.Round\${RoundCentral}.Extended.Curated.Consensus__Peak\${SelectedPeak}.fa \\
                    ./Family_\${FamilyNum}/Round_\${RoundCentral}/07_CurrentConsensi.Round\${RoundCentral}.Extended.Curated.Consensus.fa
                cp ./Family_\${FamilyNum}/Round_\${RoundCentral}/05_CurrentConsensi.Round\${RoundCentral}.Extended.Curated__Peak\${SelectedPeak}.aln.fa \\
                    ./Family_\${FamilyNum}/Round_\${RoundCentral}/07_CurrentConsensi.Round\${RoundCentral}.Extended.Curated.aln.fa

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