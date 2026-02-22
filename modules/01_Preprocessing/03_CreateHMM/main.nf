process CreateNSplitHMM {
    label "CreateHMM"

    input:
    file(AlnFile)

    output:
    path("HMM_Split*")    , emit:SplitHMM
    path("${AlnFile}.hmm"), emit:AlHMM

    script:
    """
    #Build HMM models
    hmmbuild \\
        --dna \\
        --cpu ${task.cpus} \\
        --symfrac 0 \\
        --seed 1992 \\
        --fragthresh 1.0 \\
        --wnone \\
        ${AlnFile}.hmm \\
        ${AlnFile}   

    #Split HMM models
    csplit  \\
        --suffix-format="%06d" \\
        --prefix=TEMP__HMM_Split \\
        ${AlnFile}.hmm "////+1" "{*}"

    #Remove empty files
    find . -size 0 -delete

    #Merge HMM models into bigger chunks
    Counter=1
    File=1
    for i in \$(ls TEMP__HMM_Split* | shuf);do
        if (( \$Counter % 1 == 0 ));then
            (( File++ ))
        fi
        cat \$i >> HMM_Split\$(printf "%06d" \$File)
        (( Counter++ ))
    done

    #Remove temporary files
    rm -f TEMP__HMM_Split*
    """
}