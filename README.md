# RRE – Recursive Repeat Extension

Recursive Repeat Extension (RRE) is a Nextflow DSL2 pipeline for extending and polishing repeat consensus sequences produced by tools like RepeatModeler2. It also ships an alternative HEEA workflow for a different extension strategy.

## Requirements
- Nextflow ≥ 22.03.0-edge (DSL2 enabled)
- Apptainer/Singularity or Docker images (recommended to keep dependencies consistent, we provide a Docker repository contained in dockerhub under `biofalcon/rre`, and is already included in the nextflow config file)
- Repeat libraries produced by RepeatModeler2 (`.stk` **and** `.fa` )
  
## Inputs
Required flags:
- `--Genome` — genome FASTA
- `--consensus` — consensus sequences FASTA (e.g., RepeatModeler2 output)
- `--consensusAln` — consensus alignment in Stockholm (`.stk`) alignment format
- `--workflow` — either `RRE` or `HEEA`

Optional/advanced:
- `--outDir` — output folder (default `./Results`)
- `--hmmResults` — precomputed HMMER tblout to skip running nhmmer
- `--AncientMode` — use extension strategies for older repeats
- `--hyperT` — enable hyperthreading flag passed to some tools

Key parameters (current defaults from `RRE.nf`):
- `--extension` 250 (bases added each round)
- `--maxExtensionSize` 25000
- `--maxRounds` 30 (recursive rounds)
- `--maxCentralRounds` 1 (central extension rounds)
- `--maxFamilies` 10 (per repeat family)
- `--SampleSeqs` 25 (sequences sampled per round)
- `--hmmEValue` 1e-5 (HMMER filtering)
- `--noiseThreshold` 0.2 (For alignment cleaning)
- `--prevRoundCoverage` 0.5 (Overlap that is required when sampling new sequences)
- `--percentVertical` 0.3 (vertical cleaning)
- `--polishminHits` 10 
- `--polishLengthLimitMultiplier` 1.5
- `--polishWindow` 10000
- `--polishCoverage` 50

## Quickstart
1) Clone the pipeline using Git:

```bash
git clone https://github.com/BioFalcon/RRE
```

After cloning the repository, the `netflow.config` file need to be adjusted to run in the available compute environment.

1) Run the pipeline (RRE workflow example):

```bash
nextflow run RRE.nf \
	--Genome path/to/genome.fa \
	--consensus path/to/consensi.fa \
	--consensusAln path/to/consensi.aln \
	--workflow RRE \
	--outDir ./Results \
```

Note: To run the HEEA workflow, set `--workflow HEEA` (same required inputs).

## Examples
1) Using the human genome (`HumanGenome.fa`; not included in this repository) to extend a repeat library:
```bash
nextflow run RRE.nf \
	--Genome ./HumanGenome.fa \
	--consensus ./examples/NormalRun/Consensi.fa \
	--consensusAln ./examples/NormalRun/Consensi.stk \
	--workflow RRE \
	--outDir ./Results \
```

2) Using the human genome (`HumanGenome.fa`; not included in this repository) to extend an ancient repeat (CR1_Mam):
```bash
nextflow run RRE.nf \
	--Genome ./HumanGenome.fa \
	--consensus ./examples/AncientExtension/MamCR1.Truncated.fa \
	--consensusAln ./examples/AncientExtension/MamCR1.Truncated.stk \
	--workflow RRE \
	--AncientMode \
	--SampleSeqs 150 \
	--extension 90 \
	--outDir ./Results \
```
We recommend monitoring the results of each round of extension carefully, as alignments might not reflect the reality of the extension (e.g. inclusion of other families, truncated extension). To do so, check the `Results/WorkDir/REPEAT_ID/Left/` or `RRE/Results/WorkDir/REPEAT_ID/Right` to examine each round individually. 

Sometimes the alignments are not correct fue to internal duplications , leading to incorrect merged alignments. To correct this, we recommend making a new alignment of the incorrect aligned parts using a script as follows:
```bash
mkdir MergeConsensus

cp ./Results/WorkDir/${RepeatID}/Left/Round_14/05_CurrentConsensi.Left.Round14.Extended.aln.fa ./MergeConsensus/MergedAln_Round14.aln.fa

for i in $(seq 15 1 20); do 
	RoundNum=$(echo $i | awk '{ printf("%02d\n", $1) }')
	PrevAln=$(echo $i  | awk '{ printf("%02d\n", $1 - 1) }')

	mafft --thread 32 \
		  --seed ./Results/WorkDir/${RepeatID}/Left/Round_${RoundNum}/06_CurrentConsensi.Left.Round${RoundNum}.Extended.tmp2.aln.fa \
		  --seed ./MergeConsensus/MergedAln_Round${PrevAln}.aln.fa \
		  /dev/null | \
	seqkit grep -v -r -p "DUP" > ./MergeConsensus/TEMP_Round${RoundNum}.aln.fa

	python ../../Staging/VerticalCleaning.py \
		-i ./MergeConsensus/TEMP_Round${RoundNum}.aln.fa \
		-o ./MergeConsensus/MergedAln_Round${RoundNum} \
		--perc 0.3
done
```

## Directory Structure
All files related to the run will be contained in the directory specified in the `--outDir` parameter. The structure is as follows:

<details>
<summary>For RRE</summary>

```
.
`-- Results
    `-- WorkDir
        |-- RepeatID1
        |   |-- Central - Contains the results of the first round of extensions (Module 2)
        |   |   |-- Central.FINISH.CHECK
        |   |   |-- Family_00
        |   |   |   |-- Central.GOOD.CHECK (Check to not repeat this step in case the run fails)
        |   |   |   |-- Final/
        |   |   |   |   |-- CentralPack.tar.gz 
        |   |   |   |   |-- FinalConsensi.aln.fa - Contains the alignment of the extended model 
        |   |   |   |   |-- FinalConsensi.consensus.fa - Contains the consensus sequence of the extended model
        |   |   |   |   |-- FinalConsensi.SidesCheck.Coords - Coordinates of the extended section of the repeat, the first line containes the coordinates of the 3' end, and the second the 5' end
        |   |   |   |   |-- Final.Consensus.LeftSide.Coord - Coordinates of the extended 3'section of the repeat. 
        |   |   |   |   |-- Final.Consensus.RightSide.Coord - Coordinates of the extended 5'section of the repeat. 
        |   |   |   |   |-- FinalSideCheck.aln.fa - Alignment between the consensus of the new model and the original model.
        |   |   |   |   |-- GoodCHECK - Check to skip if run is repeated.
        |   |   |   |   |-- LeftPackCentral - Directory with necesary files to extend the 3' end. 
        |   |   |   |   |   |-- FinalConsensi.aln.fa
        |   |   |   |   |   |-- FinalConsensi.consensus.fa
        |   |   |   |   |   |-- Final.Consensus.LeftSide.Coord
        |   |   |   |   |   `-- LEFT.GOOD.CHECK - Check to decide if 3' end should be extended.
        |   |   |   |   |-- LeftPack.Central.tar.gz - Compressed directory.
        |   |   |   |   |-- RightPackCentral - Directory with necesary files to extend the 3' end.
        |   |   |   |   |   |-- FinalConsensi.aln.fa
        |   |   |   |   |   |-- FinalConsensi.consensus.fa
        |   |   |   |   |   |-- Final.Consensus.RightSide.Coord
        |   |   |   |   |   `-- RIGHT.GOOD.CHECK - Check to decide if 5' end should be extended.
        |   |   |   |   `-- RightPack.Central.tar.gz - Compressed directory
        |   |   |   |-- Round_00/
        |   |   |   |   |-- 02_CurrentConsensi.Round00.Extended.bed - File containing the coordinates of individual instances.
        |   |   |   |   |-- 07_CurrentConsensi.Round00.Extended.Curated.Consensus.fa - File contanining original consensus.
        |   |   |   |   `-- GoodCHECK - Check to not generate directory again.
        |   |   |   `-- Round_01
        |   |   |       |-- 01_CurrentConsensi.Round01.Extended.bed - File containing the extended coordinates of individual isntances  ordered by bit-score
        |   |   |       |-- 02_CurrentConsensi.Round01.Extended.bed - File containing the instances that are to be used for the generation of the model  (the number is dependent of the --SampleSeqs parameter)
        |   |   |       |-- 03_CurrentConsensi.Round01.Extended.fa - Sequences to be used for the construction of the new repeat model.
        |   |   |       |-- 04_CurrentConsensi.Round01.Extended.aln - Raw alignment from the extended sequences.
        |   |   |       |-- 04_CurrentConsensi.Round01.Extended.Fams.Cluster_0.fa.aln - Alignment corresponding to Cluster 0 after the family splitting step
        |   |   |       |-- 04_CurrentConsensi.Round01.Extended.Fams.Plots.pdf - Plots corresponding to the family splitting
        |   |   |       |-- 04_CurrentConsensi.Round01.Extended.MAFFTlog
        |   |   |       |-- 05_CurrentConsensi.Round01.Extended.Curated.Consensus__Peak0.fa - Consensus of the clean alignment.
        |   |   |       |-- 05_CurrentConsensi.Round01.Extended.Curated.pdf - Plots corresponding the cleaning of the alignment.
        |   |   |       |-- 05_CurrentConsensi.Round01.Extended.Curated__Peak0.aln.fa - Clean alignment
        |   |   |       |-- 06_CurrentConsensi__Peak0.hmm - Model generated from the cleaned alignment
        |   |   |       |-- 06_CurrentConsensi__Peak0.tblout - Insntances found using the cleaned alignment
        |   |   |       |-- 07_CurrentConsensi.Round01.Extended.Curated.aln.fa - A copy of the cleaned alignment of the extended model.
        |   |   |       |-- 07_CurrentConsensi.Round01.Extended.Curated.Consensus.fa - A copy of the consensus sequence of the extended model.
        |   |   |       `-- GoodCHECK
        |   |   |-- .
        |   |   |-- .
        |   |   `-- Family_NN
        |   |-- Left
        |   |   |-- Family_00
        |   |   |   |-- GOOD.FAM.CHECK
        |   |   |   |-- Round_00
        |   |   |   |   |-- 07_CurrentConsensi.Left.Round00.Extended.Curated.aln.fa - Alignment of the extended model generated in Module 2 (Central extension).
        |   |   |   |   |-- 07_CurrentConsensi.Left.Round00.Extended.Curated.Consensus.fa - Consensus of the extended model generated in Module 2 (Central extension).
        |   |   |   |   |-- 08_CurrentConsensi.Left.Round00.Extended.Curated.hmm - HMM of the extended model from Module 2.
        |   |   |   |   |-- 08_CurrentConsensi.Left.Round00.ExtendedSide.Coord - Coordinates of the extended section of the extended model from Module 2. 
        |   |   |   |   `-- GoodCHECK - Check to not generate this directory again.
        |   |   |   |-- Round_01
        |   |   |   |   |-- 00_CurrentConsensi.Left.Round01.hmmout - nhmmer output using the previouys round model
        |   |   |   |   |-- 00_CurrentConsensi.Left.Round01.Raw.bed - Raw coordinates extracted from the nhmmer search
        |   |   |   |   |-- 01_CurrentConsensi.Left.Round01.Extended.bed - Extended coordinates of individual instances.
        |   |   |   |   |-- 02_CurrentConsensi.Left.Round01.Extended.bed - Selected extended coordinates from individual instances.
        |   |   |   |   |-- 03_CurrentConsensi.Left.Round01.Extended.fa - Genomic sequences of individual insntances to be used for repeat extension. 
        |   |   |   |   |-- 04_00_CurrentConsensi.Left.Round01.Extended.Curated.Consensus__Peak0.fa - Consensus sequence generated from cleaned alignment of the selected Family after family splitting step.
        |   |   |   |   |-- 04_00_CurrentConsensi.Left.Round01.Extended.Curated.pdf - Plots corresponding to alignment cleaning of selected Family after family splitting step.
        |   |   |   |   |-- 04_00_CurrentConsensi.Left.Round01.Extended.Curated__Peak0.aln.fa - Cleaned alignment of selected Family after family splitting step.
        |   |   |   |   |-- 04_CurrentConsensi.Left.Round01.Extended.aln - Raw alignment from the genomic sequences of the individual instances. 
        |   |   |   |   |-- 04_CurrentConsensi.Left.Round01.Extended.Fams.Cluster_0.fa.aln - Raw alignment corresponding to Cluster 0 after family splitting step.
        |   |   |   |   |-- 04_CurrentConsensi.Left.Round01.Extended.Fams.Plots.pdf - Plots corresponding to family splitting process.
        |   |   |   |   |-- 05_CurrentConsensi.Left.Round01.Extended.Curated.aln.fa - Cleaned alignment of extended sequences.
        |   |   |   |   |-- 05_CurrentConsensi.Left.Round01.Extended.Curated.Consensus.fa - Consensus sequences generated from the cleaned alignment of extended sequences.
        |   |   |   |   |-- 06_CurrentConsensi.Left.Round01.Extended.aln - Raw alignment after merging the current alignment with previous rounds' alignments.
        |   |   |   |   |-- 06_CurrentConsensi.Left.Round01.Extended.MAFFTlog 
        |   |   |   |   |-- 07_CurrentConsensi.Left.Round01.Extended.Curated.aln.fa - Cleaned merged alignment of this round.
        |   |   |   |   |-- 07_CurrentConsensi.Left.Round01.Extended.Curated.Consensus.fa - Consensus sequence generated from this round cleaned merged alignment. 
        |   |   |   |   |-- 07_CurrentConsensi.Left.Round01.Extended.Curated.pdf - Plots corresponding to the cleaning of the merged alignment.
        |   |   |   |   |-- 08_CurrentConsensi.Left.Round01.Extended.aln - Alignment between the current and previous rounds extended consensus.  
        |   |   |   |   |-- 08_CurrentConsensi.Left.Round01.Extended.Curated.hmm - HMM profile generated from the merged alignment. 
        |   |   |   |   `-- 08_CurrentConsensi.Left.Round01.ExtendedSide.Coord - Coordinates corresponding to the extended section of this round's extended model.
        |   |   |   |-- .
        |   |   |   |-- .
        |   |   |   `-- Round_NN
        |   |   |-- .
        |   |   |-- .
        |   |   |-- Family_NN
        |   |   |-- LEFT.DONE.CHECK
        |   |   `-- Staging
        |   |-- Merge
        |   |   |-- Family_00
        |   |   |   |-- FinalConsensi.aln.fa -  Merged alignment of the last round of extension.
        |   |   |   |-- FinalConsensi.aln.stk - Merged alignment of the last round of extension in stockholm format.
        |   |   |   |-- FinalConsensi.consensus.fa - Consensus generated from the alignment of the last round of extension.
        |   |   |   `-- GoodCHECK - Check not to do this step again
        |   |   `-- MERGE.GOOD.CHECK
        |   |-- Polishing
        |   |   |-- CHECKPolishing
        |   |   |-- Family_00
        |   |   |   |-- 01_ConsensusModel.hmm - HMM profile generated from the last round of extension.
        |   |   |   |-- 02_ConsensusHits.out - Results of the nhmmer search using the HMM profile from the last last round of extension.
        |   |   |   |-- 03_ConsensusHits.bed - Coordinates of individual instances.
        |   |   |   |-- 04_ConsensusHits.NoOverlaps.bed - Merged coordinates from individual instances.
        |   |   |   |-- 05_SequenceCoverage_CoveragePlot.pdf - Plots corresponding to instance selection.
        |   |   |   |-- 05_SequenceCoverage.CutCoverage.txt - Value to cut coverage at.
        |   |   |   |-- 05_SequenceCoverage.fa - Genomic sequences of the instances that will be used for the polished consensus. 
        |   |   |   |-- 05_SequenceCoverage.IndivHits.bed - 
        |   |   |   |-- 05_SequenceCoverage.Window.bed
        |   |   |   |-- 06_Polishing.Coverage.MAFFTlog
        |   |   |   |-- 06_SequenceCoverage.aln
        |   |   |   |-- 07_Polished.Curated.aln.fa
        |   |   |   |-- 07_Polished.Curated.bed
        |   |   |   |-- 07_Polished.Curated.Consensus.fa
        |   |   |   |-- 99_CoverageWarning
        |   |   |   `-- CHECKFamPolished
        |   |   `-- Staging
        |   `-- Right
        |       |-- Family_00
        |       |   |-- GOOD.FAM.CHECK
        |       |   |-- Round_00
        |       |   |-- Round_01
        |       |   |-- .
        |       |   |-- .
        |       |   `-- Round_NN
        |       |-- .
        |       |-- .
        |       |-- Family_NN
        |       |-- RIGHT.DONE.CHECK
        |       `-- Staging
        |-- .
        |-- .
        `-- RepeatIDN
```
</details>

<details>
<summary>For HEEA</summary>

```
.
`-- Results
    `-- WorkDir
        |-- RepeatID1
        |   `-- HEEA
        |       |-- HEEA.GOOD.CHECK - Check file to not re-do the extension on this repeat model. 
        |       |-- Final
        |       |   |-- Brute.tar.gz - Compressed directory with final files
        |       |   |-- FinalConsensi.aln.fa - Alignment of the last round of extension.
        |       |   |-- FinalConsensi.consensus.fa - Consensus generated from the alignment of the last round of extension
        |       |   `-- GoodCHECK - Check file to not generate this directory again.
        |       |-- Round_00
        |       |   |-- 02_CurrentConsensi.Round00.Extended.bed - Coordinates from the individual instances found in the nhmmer search.
        |       |   |-- 07_CurrentConsensi.Round00.Extended.Curated.Consensus.fa - Original consensus sequence of the repeat model. 
        |       |   `-- GoodCHECK - Check file to not generate this directory again.
        |       |-- Round_01
        |       |   |-- 01_CurrentConsensi.Round01.Extended.bed - Extended coordinates of all individual instances.
        |       |   |-- 02_CurrentConsensi.Round01.Extended.bed - Extended coordinates of individual instances that will be used for generating the repeat model.
        |       |   |-- 03_CurrentConsensi.Round01.Extended.fa - Genomic sequences of the individual instaqnces that will be used for generating the repeat model.
        |       |   |-- 04_CurrentConsensi.Round01.Extended.aln - Multiple sequence alignment of the genomic sequences.
        |       |   |-- 04_CurrentConsensi.Round01.Extended.MAFFTlog
        |       |   |-- 05_CurrentConsensi.Round01.Extended.Curated.Consensus__Peak0.fa - Consensus generated from the cleaned MSA.
        |       |   |-- 05_CurrentConsensi.Round01.Extended.Curated.pdf - Plots corresponding to MSA cleaning.
        |       |   |-- 05_CurrentConsensi.Round01.Extended.Curated__Peak0.aln.fa - Cleaned MSA.
        |       |   |-- 07_CurrentConsensi.Round01.Extended.Curated.aln.fa - A copy of the cleaned MSA.
        |       |   |-- 07_CurrentConsensi.Round01.Extended.Curated.Consensus.fa -  A copy of the consensus generated from the cleaned MSA.
        |       |   |-- 08_CurrentConsensi.Round01.Consensi.aln.fa - Alignment between the current and previous consensus sequence.
        |       |   |-- 08_CurrentConsensi.Round01.Consensi.coord.bed - Coordinates corresponding to the regions that have been extended in this round. The first line corresponds to the 3' end, and the second line to the 5' end.
        |       |   `-- GoodCHECK - Check file to not generate this directory again.
        |       |-- .
        |       |-- .
        |       |-- Round_NN
        |       `-- Staging
        |-- .
        |-- .
        `-- RepeatIDN
```
</details>

## Outputs
- `Results/WorkDir/` — per-repeat working directories containing intermediate HMMER hits, alignments, curated consensi, and logs.
- `Results/MergingLibrary` — final polished consensus sequences and checks from the Polishing and CD-HIT steps. Results are provided in `.stk` and `.fa` formats

## Notes & tips
- If you already have HMMER search results, provide `--hmmResults` to skip `HMMER_Run`. This saves time in case you have already ran RRE before.
- If extending a custom repeat library, make sure the names in the `.fasta` and `.stk` files match. The ID in the `.stk` file is contained in the `#=GF ID` field. In the fasta file, make sure the ID is followed by a `#` character, otherwise the pipeline won't work (**IMPORTANT**)
- If extending a model from DFAM (or other libraries that may contain very gappy alignments), we recommend running nhmmer using the HMM model provided by DFAM separately, and providing the results to RRE using the `--hmmResults` flag.
- Modify the provided nextflow.config file to take advantage of your HPC environment. We provide some values for memory, CPUs, and walltime, but this can be modified according to the available resources. While we use `slurm` as template, other workflow managers such as `pbs` and `sge` are compatible with Nextflow.
- Adjust to a higher `--extension` value to do less rounds, but increasing it might lead to chimeric families (e.g. inclussion of Alu families in other repeat families)
- RRE was built with "resumability" in mind. In case Nextflow is unable to resume a run through the `-resume` flag, as longs as the `--outDir` remains the same RRE will not overwrite any of the already generated results and will resume at the latest checkpoint.
- RRE is agnostic to the input genome. We observed long runtimes when segmental duplications are present. Usage of tools like [BISER](https://github.com/0xTCG/biser) can help detect this duplications to hard-mask them before starting a run.
- We like to use [AliView](https://ormbunkar.se/aliview/) to visualize our MSAs.
  
## Citation
Recursive Repeat Extender (RRE): A recursive approach to automatically extend repeat element models
Francisco Falcon, Elly M. Tanaka, Diego Rodriguez-Terrones
bioRxiv 2026.04.14.718546; doi: https://doi.org/10.64898/2026.04.14.718546
