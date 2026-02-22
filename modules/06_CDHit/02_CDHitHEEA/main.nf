process CDHitHEEA{
    label "CDHit"

    input:
    path("./InputFiles/FINALCHECK*")
    path(outDir)
    path(ConsensusFasta)
    path(ConsensusAln)

    script:
    """
    #Save workdir into variable
    oldWorkDir=\$(pwd)

    #Move to new WorkDir 
    mkdir -p ${outDir}/MergingLibrary
    cd ${outDir}/MergingLibrary

    #Backup new path 
    NewPath=\$(pwd)

    #Copy files to new path
    cp \${oldWorkDir}/${ConsensusFasta} ./
    cp \${oldWorkDir}/${ConsensusAln} ./

    #Concatenate all HEEA extended consensi
    cat ../WorkDir/*/HEEA/Final/FinalConsensi.consensus.fa > Extended.fa

    #Extract sequences that were not extended
    ##Make fai index of extended sequences
    seqkit faidx Extended.fa

    ##Extract non-extended sequences
    seqkit grep -v -r -f <(cut -f1 Extended.fa.fai | sed 's/#.*/#/') ${ConsensusFasta} > NonExtended.fa

    #Concatenate extended and non-extended sequences
    cat Extended.fa NonExtended.fa > RedundantConsensi.fa

    #Run CD-Hit
    ##Calculate available memory
    AvailMem=\$( echo ${task.memory.toMega()} | awk '{print int(\$1*0.90)}' )

    ##CD-Hit
    cd-hit-est \\
    -i RedundantConsensi.fa \\
    -o RedundantConsensi.fa.merged \\
    -G 0 -b 500 \\
    -M \${AvailMem} \\
    -T ${task.cpus} \\
    -c 0.8 \
    -aS 0.8

    #Extract stk sequences
    > RedundantConsensi.fa.merged.stk

    ##For extended sequences
    for i in \$(cut -f1 Extended.fa.fai| sed 's/#.*//');do
        echo "# STOCKHOLM 1.0" >> RedundantConsensi.fa.merged.stk
        echo -e "#=GF ID\\t\${i}" >> RedundantConsensi.fa.merged.stk
        cat ../WorkDir/\${i}/HEEA/Final/FinalConsensi.aln.fa | \\
            seqkit replace -p "-" -s -r "." | \\
            seqkit fx2tab | \\
            cut -f1,2 | \\
            sed 's/)\\t/\\t/;s/(/_/' \\
        >> RedundantConsensi.fa.merged.stk
        echo "//" >> RedundantConsensi.fa.merged.stk
    done

    ##For non-extended sequences
    
    seqkit faidx RedundantConsensi.fa.merged
    for i in \$(grep -v -f <(cut -f1 Extended.fa.fai ) RedundantConsensi.fa.merged.fai | cut -f1);do
        ID=\$(echo \$i | sed 's/#.*//')
        echo "# STOCKHOLM 1.0" >> RedundantConsensi.fa.merged.stk
        sed -n "/\${ID}\$/,/\\/\\//p" ${ConsensusAln} >> RedundantConsensi.fa.merged.stk
    done

    #Annotate CD-Hit output
    RepeatClassifier -consensi ./RedundantConsensi.fa.merged \\
                     -threads ${task.cpus} \\
                     -stockholm ./RedundantConsensi.fa.merged.stk

    #Copy results to oldDir
    ln -s \${PWD}/ \${oldWorkDir}/

    #Go back to old workdir
    cd \${oldWorkDir}
    """
}