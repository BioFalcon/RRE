process FAIndex{
    label 'preprosessing_low'

    input:
    file(Genome)
    
    output:
    file("${Genome}.fai")

    script:
    """
    samtools faidx ${Genome}
    """
}