process SeedPrep {
    label "CreateHMM"

    input:
    file(AlnFile)
    file(SeedScript)

    output:
    path("Seeds.stk")         , emit:Stk
    path("Seeds.consensus.fa"), emit:Fasta

    script:
    RFMode = params.rfAll ? "--rfAll" : params.rfFromOccupancy ? "--rfFromOccupancy --rfMinOcc ${params.rfMinOcc}" : ""
    """
    #Match columns of each seed family (its #=GC RF line, or --rfAll / --rfFromOccupancy) and its class
    python3 ${SeedScript} \\
        --input ${AlnFile} \\
        --output Seeds \\
        ${RFMode}

    #Consensus of each family as AncientMode builds its starting model: from the match columns (--hand)
    hmmbuild \\
        --hand \\
        --dna \\
        --informat stockholm \\
        --seed 1992 \\
        --cpu ${task.cpus} \\
        --fragthresh 1.0 \\
        --maxinsertlen 10 \\
        --eentexp \\
        Seeds.hmm \\
        Seeds.stk > /dev/null

    #Name them ID#Class, as the rest of the pipeline expects
    hmmemit -c Seeds.hmm | \\
        awk 'NR==FNR{Class[\$1]=\$2; next} /^>/{Name=substr(\$1,2); sub(/-consensus\$/,"",Name); print ">"Name"#"Class[Name]; next}{print toupper(\$0)}' \\
        Seeds.classes.tsv - > Seeds.consensus.fa
    """
}
