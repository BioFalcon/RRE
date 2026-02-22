process Polishing {
    label 'Polishing'

    input:
    tuple val(RepeatID), path(MergeCHECK)
    path(outDir)
    path(CoverageScript)
    path(CurateScript)
    path(GenomeFasta)
    path(HMMDB_Genome)

    output:
    tuple val(RepeatID), path("CHECKPolishing")

    script:
    """
    #Load ID of the repeat into a variable
    Repeat=\$(echo ${RepeatID}| sed 's/#.*//')

    #Save workdir into variable
    oldWorkDir=\$(pwd)

    #Move to new WorkDir 
    cd ${outDir}/WorkDir/\${Repeat}

    #Create directory for repeat
    mkdir -p ./Polishing
    cd ./Polishing

    #Backup new path 
    NewPath=\$(pwd)

    #Make staging directory
    mkdir -p ./Staging

    #Copy files into work dir
    ln -fs \${oldWorkDir}/* ./Staging/

    #Check if the polishing was already done
    if [[ ! -f CHECKPolishing ]]; then
        for CurrFamily in \$( ls -d ../Merge/Family*| sed 's/.*_//' );do
            #Obtain consensus length
            consLen=\$(seqkit fx2tab --length --name  ../Merge/Family_\${CurrFamily}/FinalConsensi.consensus.fa | cut -f2)

            #Make the directory 
            if [[ ! -d ./Family_\${CurrFamily}/ ]];then
                mkdir -p ./Family_\${CurrFamily}/
            fi

            #Check if family is already polished and contains sequence
            if [[ ! -f ./Family_\${CurrFamily}/CHECKFamPolished ]];then
                if [[ \${consLen} -gt 10 ]]; then

                    if [[ ! -f ./Family_\${CurrFamily}/01_HMMERCheck ]];then
                        #Remove old files
                        if [[ -d ./Family_\${CurrFamily}/ ]];then
                            rm -r ./Family_\${CurrFamily}/
                            mkdir -p ./Family_\${CurrFamily}/
                        fi

                        #Search the consensus using HMMER
                        ##Create HMM profile
                        hmmbuild \
                            --symfrac 0 \
                            --dna \
                            --informat afa \
                            --seed 1992 \
                            --fragthresh 1.0 \
                            --wnone \
                            -n \${Repeat}__\${CurrFamily} \
                            ./Family_\${CurrFamily}/01_ConsensusModel.hmm \
                            ../Merge/Family_\${CurrFamily}/FinalConsensi.aln.fa

                        ##Search the consensus
                        nhmmer \
                            --cpu ${task.cpus} \
                            --tblout ./Family_\${CurrFamily}/02_ConsensusHits.out \
                            --notextw \
                            --noali \
                            --seed 1992 \
                            -E 1e-5 \
                            ./Family_\${CurrFamily}/01_ConsensusModel.hmm \
                            ./Staging/HMMDB_Genome
                    else
                        echo "HMMER already run, skipping to next step"
                    fi

                    if [[ \$? -eq 0 ]];then
                        > ./Family_\${CurrFamily}/01_HMMERCheck
                    fi

                    #Check if there are enough hits
                    if [[ \$(sed '/^#/d' ./Family_\${CurrFamily}/02_ConsensusHits.out | wc -l ) -gt ${params.polishminHits} ]];then 
                        #Convert HMM output to bed
                        awk '{OFS="\\t"}
                            {if (\$12 == "+"){
                                print \$1,\$9,\$10,\$5"-"\$6,\$14,\$12
                            }else{
                                print \$1,\$10,\$9,\$5"-"\$6,\$14,\$12
                            } 
                            }' ./Family_\${CurrFamily}/02_ConsensusHits.out > \\
                        ./Family_\${CurrFamily}/03_ConsensusHits.bed

                        #Perform BedMap
                        ##On the possitive strand
                        awk '\$6=="+"{print}' ./Family_\${CurrFamily}/03_ConsensusHits.bed | \\
                        sort -k1,1 -k2,2n | \\
                        bedmap --count --echo-map-range --echo-map-id --echo-map-score  --echo-map-size --fraction-either 1 --delim '\\t' - | \\
                        cut -f2- | \\
                        sort -k1,1 -k2,2n | \\
                        uniq | \\
                        awk '{OFS="\\t";FS="\\t"}{ \\
                                split(\$5, scores, ";"); \\
                                split(\$4, order, ";"); \\
                                max_score = -999999999; \\
                                max_index = 0; \\
                                for (i = 1; i <= length(scores); i++) { \\
                                    if (scores[i] > max_score) { \\
                                        max_score = scores[i]; \\
                                        max_index = order[i]; \\
                                    } \\
                                } \\
                                \$5 = max_score; \
                                \$4 = max_index; \
                                print \$1, \$2, \$3, \$4, \$5, "+";}' | \
                        sort -k1,1 -k2,2n | \\
                        uniq | \\
                        sed "s/\$/\\t/" > ./Family_\${CurrFamily}/04_ConsensusHits.NoOverlaps.bed

                        ##On the negative strand
                        awk '\$6=="-"{print}' ./Family_\${CurrFamily}/03_ConsensusHits.bed  | \\
                        sort -k1,1 -k2,2n | \\
                        bedmap --count --echo-map-range --echo-map-id --echo-map-score  --echo-map-size --fraction-either 1 --delim '\\t' - | \\
                        cut -f2- | \\
                        sort -k1,1 -k2,2n | \\
                        uniq | \\
                        awk '{OFS="\\t";FS="\\t"}{ \\
                            split(\$5, scores, ";"); \\
                            split(\$4, order, ";"); \\
                            max_score = -999999999; \\
                            max_index = 0; \\
                            for (i = 1; i <= length(scores); i++) { \\
                                if (scores[i] > max_score) { \\
                                    max_score = scores[i]; \\
                                    max_index = order[i]; \\
                                } \\
                            } \\
                            \$5 = max_score; \\
                            \$4 = max_index; \\
                            print \$1, \$2, \$3, \$4, \$5, "-";}'  | \\
                        sort -k1,1 -k2,2n | \\
                        uniq | \\
                        sed "s/\$/\\t/" >> ./Family_\${CurrFamily}/04_ConsensusHits.NoOverlaps.bed

                        #Get the sequences for the new alignment
                        ./Staging/${CoverageScript} \\
                        --input ./Family_\${CurrFamily}/04_ConsensusHits.NoOverlaps.bed \\
                        --output ./Family_\${CurrFamily}/05_SequenceCoverage \\
                        --genome ./Staging/${GenomeFasta} \\
                        --size \${consLen} \\
                        --target ${params.polishCoverage} \\
                        --window ${params.polishWindow}

                        #Define if sensitive alignment can be done
                        MafftPars="--localpair"

                        #Check if coverage went overboard
                        ##Calculate median 
                        MedianCov=\$(seqkit stats -T ./Family_\${CurrFamily}/05_SequenceCoverage.fa | \\
                                    cut -f7 | \\
                                    sed '1d' | \\
                                    awk '{print int(\$1)}' )

                        if [[ \$MedianCov -gt \$( echo \${consLen} | awk '{print int(\$1*${params.polishLengthLimitMultiplier})}' ) ]];then
                            #Just pass on the previous consensus, make check and note
                            cp ../Merge/Family_\${CurrFamily}/FinalConsensi.consensus.fa ./Family_\${CurrFamily}/07_Polished.Curated.Consensus.fa
                            cp ../Merge/Family_\${CurrFamily}/FinalConsensi.aln.fa ./Family_\${CurrFamily}/07_Polished.Curated.aln.fa
                            touch ./Family_\${CurrFamily}/07_Polished.Curated.bed
                            echo "Coverage went overboard, passing on previous consensus" > ./Family_\${CurrFamily}/99_CoverageWarning
                        else
                            #Perform alignment
                            mafft \\
                            \${MafftPars} \\
                            --maxiterate 1000 \\
                            --thread ${task.cpus} \\
                            ./Family_\${CurrFamily}/05_SequenceCoverage.fa \\
                            >  ./Family_\${CurrFamily}/06_SequenceCoverage.aln \\
                            2> ./Family_\${CurrFamily}/06_Polishing.Coverage.MAFFTlog
                            
                            #Setup variables for curation
                            ThesLen=\$(cat ./Family_\${CurrFamily}/05_SequenceCoverage.CutCoverage.txt)

                            #Curate alignment
                            ./Staging/${CurateScript} \\
                            --input ./Family_\${CurrFamily}/06_SequenceCoverage.aln \\
                            --output ./Family_\${CurrFamily}/07_Polished.Curated \\
                            --SeqID \${Repeat}#__Family\${CurrFamily}__Polished \\
                            --zeroOnly

                            #move some files
                            mv ./Family_\${CurrFamily}/07_Polished.Curated__Peak0.aln.fa ./Family_\${CurrFamily}/07_Polished.Curated.aln.fa
                            mv ./Family_\${CurrFamily}/07_Polished.Curated.Consensus__Peak0.fa ./Family_\${CurrFamily}/07_Polished.Curated.Consensus.fa

                            #If new consensus went overboard, pass on the previous consensus
                            if [[ \$(seqkit stats -T ./Family_\${CurrFamily}/07_Polished.Curated.Consensus.fa | cut -f5 | sed '1d') -gt \$( echo \${consLen} | awk '{print int(\$1*${params.polishLengthLimitMultiplier})}' ) ]];then
                                cp ../Merge/Family_\${CurrFamily}/FinalConsensi.consensus.fa ./Family_\${CurrFamily}/07_Polished.Curated.Consensus.fa
                                cp ../Merge/Family_\${CurrFamily}/FinalConsensi.aln.fa ./Family_\${CurrFamily}/07_Polished.Curated.aln.fa
                                touch ./Family_\${CurrFamily}/07_Polished.Curated.bed
                                echo "Coverage went overboard, passing on previous consensus" > ./Family_\${CurrFamily}/99_CoverageWarning
                            fi
                            
                            #If consensus is empty, add a base so it doesnt fail
                            if [[ \$(seqkit stats -T ./Family_\${CurrFamily}/07_Polished.Curated.Consensus.fa | cut -f5 | sed '1d') -eq 0 ]];then
                                echo "A" >> ./Family_\${CurrFamily}/07_Polished.Curated.Consensus.fa
                            fi
                            cp ./Family_\${CurrFamily}/05_SequenceCoverage.IndivHits.bed ./Family_\${CurrFamily}/07_Polished.Curated.bed
                        fi

                    else
                        cp ../Merge/Family_\${CurrFamily}/FinalConsensi.consensus.fa ./Family_\${CurrFamily}/07_Polished.Curated.Consensus.fa
                        cp ../Merge/Family_\${CurrFamily}/FinalConsensi.aln.fa ./Family_\${CurrFamily}/07_Polished.Curated.aln.fa
                        touch ./Family_\${CurrFamily}/07_Polished.Curated.bed
                        echo "Coverage went overboard, passing on previous consensus" > ./Family_\${CurrFamily}/99_CoverageWarning
                    fi
                    
                else
                    #If consensus is empty, create empty files
                    cp ../Merge/Family_\${CurrFamily}/FinalConsensi.consensus.fa ./Family_\${CurrFamily}/07_Polished.Curated.Consensus.fa
                    echo A >> ./Family_\${CurrFamily}/05_Polishing.Curated.Consensus.fa 
                    cp ../Merge/Family_\${CurrFamily}/FinalConsensi.consensus.fa ./Family_\${CurrFamily}/07_Polished.Curated.aln.fa
                    touch ./Family_\${CurrFamily}/07_Polished.Curated.bed
                    #Make waring file
                    echo "Consensus is empty, passing on previous consensus" > ./Family_\${CurrFamily}/99_CoverageWarning
                fi
            fi

            #Make a check file
            touch ./Family_\${CurrFamily}/CHECKFamPolished
        done

        #Make a check file
        touch CHECKPolishing
    fi 

    #Copy results to oldDir
    ln -s \${PWD}/CHECKPolishing \${oldWorkDir}/

    #Go back to old workdir
    cd \${oldWorkDir}
    """
}