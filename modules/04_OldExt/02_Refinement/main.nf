process RRE_Refinement {
    label "RRE_Refinement"

    input:
    tuple val(RepeatID), path(MergeCHECK)
    path(HMMRDB)
    path(GenomeIndex)
    path(Genome)
    path(outDir)
    path(RefineScript)

    output:
    tuple val(RepeatID), path("REFINE.CHECK")

    script:
    """
    #Genome-wide refinement of the merged model (AncientMode). It is its own step because nhmmer with the
    #full-length model against the whole genome needs much more memory than the extension rounds. If it is
    #killed (memory or time limit), the step fails and Nextflow retries it with more resources; the extension
    #is not repeated.
    #The final model of the repeat goes to FinalModel/, next to Merge/: the refined model, or the merged one
    #(with NOT_REFINED.txt saying why) if the refinement is skipped or fails.

    #Save nf directory into a variable
    oldWorkDir=\$(pwd)

    #Move to new WorkDir
    cd ${outDir}/WorkDir/${RepeatID}
    mkdir -p ./FinalModel
    cd ./FinalModel

    #Done in a previous run: skipped, unless the merged model changed since (the extension was redone)
    if [ -f ./REFINE.CHECK ] && [ ../Merge/Merge.aln.fa -nt ./REFINE.CHECK ];then
        echo "The merged model is newer than the final model: refining again"
        rm ./REFINE.CHECK
    fi

    if [ ! -f ./REFINE.CHECK ];then
        #Start clean: a previous attempt may have been killed halfway
        rm -rf ./${RepeatID}.* ./Iterations ./NOT_REFINED.txt
        Reason=""

        if [ ${params.refineIterations} -eq 0 ];then
            Reason="Refinement skipped (--refineIterations 0)"
        elif [ ! -f ../Merge/Merge.aln.fa ];then
            Reason="No merged model to refine"
        else
            echo "Refining the final model genome-wide (${params.refineIterations} iterations at most)"
            Status=0
            bash \${oldWorkDir}/${RefineScript} \\
                -a ../Merge/Merge.aln.fa \\
                -d \${oldWorkDir}/${HMMRDB} \\
                -g \${oldWorkDir}/${Genome} \\
                -o ./${RepeatID} \\
                -W ./Iterations \\
                -n ${RepeatID} \\
                -t ${task.cpus} \\
                -I ${params.refineIterations} \\
                -Q ${params.rfMinOcc} \\
                -D ${params.maxIdentity} || Status=\$?

            #Killed by a signal (out of memory, time limit): fail, so that Nextflow retries with more resources
            if [ \${Status} -gt 128 ];then
                echo "Refinement killed (exit \${Status})"
                exit \${Status}
            #Not run at all (script or program not found): an error of the pipeline, stop
            elif [ \${Status} -eq 126 ] || [ \${Status} -eq 127 ];then
                echo "Refinement could not run (exit \${Status})"
                exit \${Status}
            #Any other failure (e.g. too few hits): keep the merged model
            elif [ \${Status} -ne 0 ];then
                Reason="Refinement failed (exit \${Status}), see the log of the RRE_Refinement step"
            fi
        fi

        #Not refined: the merged model is the final one
        if [ -n "\${Reason}" ];then
            echo "\${Reason}; the final model is the merged model" | tee ./NOT_REFINED.txt
            rm -rf ./${RepeatID}.*
            if [ -f ../Merge/Merge.aln.fa ];then
                sed "1s/.*/>${RepeatID}/" ../Merge/Merge.Consensus.fa > ./${RepeatID}.Consensus.fa
                cp ../Merge/Merge.aln.fa ./${RepeatID}.aln.fa
            fi
        fi

        #Make check
        > ./REFINE.CHECK
    fi

    #Link results to oldDir
    ln -s \$( realpath ./REFINE.CHECK ) \${oldWorkDir}/
    cd \${oldWorkDir}
    """
}
