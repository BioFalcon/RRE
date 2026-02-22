process HMMDBBuild {
    label "HMMDBBuild"

    input:
    file(Genome)

    output:
    file("HMMDB_Genome")

    script:
    """
    #Build HMM database
    makehmmerdb \\
        --informat fasta \\
        ${Genome} \\
        HMMDB_Genome
    """
}