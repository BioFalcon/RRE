process CDHit{
    label "CDHit"

    input:
    path("./InputFiles/PolishInput*")
    path(outDir)
    path(ConsensusFasta)
    path(ConsensusAln)

    output:

    script:
    """
    #Save workdir into variable
    oldWorkDir=\$(pwd)

    #Move to new WorkDir 
    mkdir -p ${outDir}/MergingLibrary
    cd ${outDir}/MergingLibrary

    #Backup new path 
    NewPath=\$(pwd)

    #Copy consensus fasta and aln to new workdir
    cp \${oldWorkDir}/${ConsensusFasta} \${oldWorkDir}/${ConsensusAln} ./
    
    #Make table with substitutions
    ls ../WorkDir/*/Polishing/Family_*/07_Polished.Curated.Consensus.fa | \\
        sed 's/.*WorkDir\\///;s/\\/Polishing\\//\\t/;s/\\/.*//'| awk '{OFS="\\t"}{printf("%s\\tSeq_%05d\\n",\$0,NR)}' \\
    > ExtendedIDs.table

    #Make fasta file with new IDs
    ID=00001
    > Extended.fa
    for i in \$(ls ../WorkDir/*/Polishing/Family_*/07_Polished.Curated.Consensus.fa); do
        seqkit replace -p "(.+)" -r "Seq_\${ID}" \$i >> Extended.fa
        ID=\$(echo \$ID | awk '{printf("%05d", \$1+1)}')
    done

    #Extract sequences that were not extended
    seqkit grep -v -r -f <(cut -f1 ExtendedIDs.table| sed 's/\$/#/') ${ConsensusFasta} > NonExtended.fa

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
    -c 0.8 \\
    -aS 0.8

    #Extract stk files
    ##Make index
    seqkit faidx RedundantConsensi.fa.merged
    seqkit faidx Extended.fa

    > RedundantConsensi.fa.merged.stk
    ##For extended sequences
    for i in \$(grep Seq_ RedundantConsensi.fa.merged.fai| cut -f1);do
        Family=\$(grep \$i ExtendedIDs.table| cut -f2)
        ID=\$(grep \$i ExtendedIDs.table| cut -f1)
        echo "# STOCKHOLM 1.0" >> RedundantConsensi.fa.merged.stk
        echo -e "#=GF ID\\t\${i}" >> RedundantConsensi.fa.merged.stk
        cat ../WorkDir/\${ID}/Polishing/\${Family}/07_Polished.Curated.aln.fa | \\
            seqkit replace -p "-" -s -r "." | \\
            seqkit fx2tab | \\
            cut -f1,2 | \\
            sed 's/__Window.*\\t/\\t/;s/)\\t/\\t/;s/(/_/' \\
        >> RedundantConsensi.fa.merged.stk
        echo "//" >> RedundantConsensi.fa.merged.stk
    done

    ##For non-extended sequences
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
    ln -s \${PWD}/RedundantConsensi.fa.merged.stk \${oldWorkDir}/
    ln -s \${PWD}/RedundantConsensi.fa.merged \${oldWorkDir}/
    ln -s \${PWD}/RedundantConsensi.fa.merged.classified \${oldWorkDir}/
    ln -s \${PWD}/RedundantConsensi.fa.merged-classified.stk \${oldWorkDir}/

    #Go back to old workdir
    cd \${oldWorkDir}
    """

}