process RRE_RecursiveExtension {
    label "RRE_Extension"

    input:
    tuple val(RepeatID),file(CentralCHECK)
    file(HMMRDB)
    file(GenomeIndex)
    file(outDir)
    file(CurateScript)
    file(OutlierScript)
    file(Genome)
    file(RREScript)
    file(ExtractCoordScript)
    file(VerticalScript)

    output:
    tuple val(RepeatID), path("MERGE.GOOD.CHECK")

    script:
    """
    bash ${RREScript} \\
        -i ${RepeatID} \\
        -t ${task.cpus} \\
        -e ${params.extension} \\
        -g ${Genome} \\
        -h ${HMMRDB} \\
        -x ${GenomeIndex} \\
        -o ${outDir} \\
        -C ${CurateScript} \\
        -O ${OutlierScript} \\
        -s ${ExtractCoordScript} \\
        -m ${params.maxExtensionSize} \\
        -p ${params.prevRoundCoverage} \\
        -S ${params.SampleSeqs} \\
        -R ${params.maxRounds} \\
        -T ${params.noiseThreshold} \\
        -P ${params.percentVertical} \\
        -V ${VerticalScript}
    """

}