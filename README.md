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
