process HMMER_Run {
    label "HMMER_Run"

    input:
    file(HMMDB)
    file(HMM_Models)

    output:
    path("${HMM_Models}.Filtered.tblout")

    script:
    """
    #Run nhmmer
    nhmmer \\
        --cpu ${task.cpus} \\
        --tblout ${HMM_Models}.tblout \\
        --notextw \\
        --noali \\
        --seed 1992 \\
        ${HMM_Models} \\
        HMMDB_Genome

    #Remove headers and change whitespace to tab
    sed -i '/^#/d;s/\\s\\{1,\\}/\\t/g' ${HMM_Models}.tblout

    #Filter Results
    awk ' \$13 <= ${params.hmmEValue} {print}' ${HMM_Models}.tblout > ${HMM_Models}.Filtered.tblout
    """
}