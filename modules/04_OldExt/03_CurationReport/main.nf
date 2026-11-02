process RRE_CurationReport {
    label "preprosessing_low"

    input:
    tuple val(RepeatID), path(RefineCHECK)
    path(outDir)
    path(ReportScript)

    output:
    tuple val(RepeatID), path("REPORT.DONE")

    script:
    """
    #Report of the alignment cleaning in the extension rounds (AncientMode): for every round and failed attempt,
    #the distribution of the column scores of the raw alignment, the valley selected, the cap and the threshold used
    #(from the Curation.json files the rounds write); then the final model from the alignment of its copies (depth,
    #agreement with the consensus, deletions, call margin, where copies start and end) and the same metrics for the
    #merged alignment and each refinement iteration, overlaid. Quick, so it is made again whenever the step runs.
    #A failed report does not stop the run.
    mkdir -p ${outDir}/WorkDir/${RepeatID}/Report
    python3 ${ReportScript} \\
        --workdir ${outDir}/WorkDir/${RepeatID} \\
        --out ${outDir}/WorkDir/${RepeatID}/Report/${RepeatID} \\
        || echo "The curation report of ${RepeatID} could not be made (exit \$?)"
    > ./REPORT.DONE
    """
}
