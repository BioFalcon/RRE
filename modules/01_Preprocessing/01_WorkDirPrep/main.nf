process WorkDirPrep{
    label 'local'

    input:
    path(outDir)

    script:
    """
    mkdir -p ${outDir}/WorkDir
    """
}
