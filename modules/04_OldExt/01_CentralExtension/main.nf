process RRE_CentralExtension{
    label "RRE_Extension"

    input:
    file(Consensus)
    file(HMMROut)
    file(outDir)
    file(SeedAln)

    output:
    tuple env(Repeat),file("Central.FINISH.CHECK"), emit:AllSides

    script:
    """
    #AncientMode has no central extension: the seed is the starting model and both sides extend from it.
    #This step builds that model from the seed's match columns and hands it, with the family's genomic
    #hits, to the recursive extension of each side.

    #Load ID of the repeat into a variable
    Repeat=\$( cat ${Consensus} )

    #Save nf directory into a variable
    oldWorkDir=\$(pwd)

    #Move to new WorkDir
    cd ${outDir}/WorkDir
    mkdir -p \${Repeat}/Central
    cd \${Repeat}/Central

    ExtensionSize=${params.extension}
    Round=./Family_00/Round_01
    Final=./Family_00/Final

    if [ ! -f \${Final}/GoodCHECK ];then
        rm -rf ./Family_00
        mkdir -p \${Round} \${Final}

        #Genomic hits of the seed in the first genome search, best first
        awk -v Repeat="\${Repeat}" '{OFS="\\t"} \$3==Repeat{print}' \${oldWorkDir}/${HMMROut} | \\
            sort -k14,14nr | \\
            cut -f1,9,10,12,14 | \\
            awk '{OFS="\\t"}{ if (\$4 == "+"){print \$1,\$2,\$3,".",\$5,\$4}else{print \$1,\$3,\$2,".",\$5,\$4}}' | \\
            sort -k5,5nr > \${Round}/02_CurrentConsensi.Round01.Extended.bed

        #This family's block of the prepared seed alignment (SeedPrep: match columns set, one block per family)
        awk -v ID="\${Repeat}" '{Block=Block \$0 "\\n"} /^\\/\\//{if (index(Block, "\\n#=GF ID " ID "\\n")) printf "%s", Block; Block=""}' \\
            \${oldWorkDir}/${SeedAln} > \${Round}/Seed.stk
        [ -s \${Round}/Seed.stk ] || { echo "ERROR: \${Repeat} not found in the prepared seed alignment"; exit 1; }

        #The starting model, from the seed's match columns (#=GC RF), as Dfam builds its models
        hmmbuild --hand --dna --informat stockholm --seed 1992 --cpu ${task.cpus} --fragthresh 1.0 --maxinsertlen 10 --eentexp -n \${Repeat} \\
            \${Round}/07_CurrentConsensi.Round01.Extended.Curated.hmm \${Round}/Seed.stk > /dev/null || \\
            { echo "ERROR: could not build the model of the seed of \${Repeat}"; exit 1; }
        hmmemit -c \${Round}/07_CurrentConsensi.Round01.Extended.Curated.hmm | \\
            awk '/^>/{print; next}{print toupper(\$0)}' \\
            > \${Round}/07_CurrentConsensi.Round01.Extended.Curated.Consensus.fa

        #The seed's match columns only, one column per model position
        python3 - \${Round}/Seed.stk > \${Round}/07_CurrentConsensi.Round01.Extended.Curated.aln.fa <<'EOF'
import sys
Seqs, RF, Order = {}, "", []
for line in open(sys.argv[1]):
    line = line.rstrip("\\n")
    if not line or line.startswith("//") or line.startswith("# STOCKHOLM"):
        continue
    if line.startswith("#=GC RF"):
        RF += line.split()[-1]
    elif not line.startswith("#"):
        Name, Seq = line.split()[:2]
        if Name not in Seqs:
            Seqs[Name] = []; Order.append(Name)
        Seqs[Name].append(Seq)
Keep = [j for j, c in enumerate(RF) if c not in ".-~_"]
for n in Order:
    Row = "".join("".join(Seqs[n])[j] for j in Keep).upper().replace(".", "-")
    if Row.strip("-"):
        print(">" + n + "\\n" + Row)
EOF
        ModelLen=\$(grep "^LENG" \${Round}/07_CurrentConsensi.Round01.Extended.Curated.hmm | awk '{print \$2}')
        echo "Seed used as the starting model: \${ModelLen} positions, \$(grep -c '>' \${Round}/07_CurrentConsensi.Round01.Extended.Curated.aln.fa) sequences, \$(wc -l < \${Round}/02_CurrentConsensi.Round01.Extended.bed) genomic hits"

        #Hand the model to both sides: each extends from it and stops by itself if it finds no extension
        cp \${Round}/07_CurrentConsensi.Round01.Extended.Curated.aln.fa \${Final}/FinalConsensi.aln.fa
        cp \${Round}/07_CurrentConsensi.Round01.Extended.Curated.Consensus.fa \${Final}/FinalConsensi.consensus.fa
        cp \${Round}/07_CurrentConsensi.Round01.Extended.Curated.hmm \${Final}/FinalConsensi.hmm
        > \${Final}/FinalConsensi.SeedModel
        cp \${Round}/02_CurrentConsensi.Round01.Extended.bed \${Final}/FinalConsensi.Extension.Coord.bed
        #Search windows: the outer ExtensionSize positions of each end of the model
        echo -e "1\\t\$(( ExtensionSize < ModelLen ? ExtensionSize : ModelLen ))" > \${Final}/Final.Consensus.LeftSide.Coord
        echo -e "\$(( ModelLen > ExtensionSize ? ModelLen - ExtensionSize + 1 : 1 ))\\t\${ModelLen}" > \${Final}/Final.Consensus.RightSide.Coord

        > \${Final}/GoodCHECK
    fi

    #Make Check for finished step
    > ./Central.FINISH.CHECK

    #Link results to oldDir
    ln -s \$( realpath ./Central.FINISH.CHECK ) \${oldWorkDir}
    cd \${oldWorkDir}
    """
}
