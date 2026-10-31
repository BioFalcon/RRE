# Changelog

## AncientModeOptimization (unreleased)

Branched from `dev_mafft_profile` (42e0993). Every change below affects only the
`--AncientMode` extension unless stated otherwise; the normal RRE and HEEA workflows
are unchanged. New command line options are listed in the README under
"AncientMode options".

### Starting model
- **The seed alignment is the only input in AncientMode; `--consensus` is ignored.** A new
  first step (`SeedPrep`, `utils/SeedPrep.py`) reads each family's ID (`#=GF ID`), class
  (`#=GF TP`, `Unknown` if absent) and match columns, and derives the consensus from the model
  built on those columns (`hmmbuild --hand`, `hmmemit -c`), so it always matches the starting
  model position by position.
- **Match columns must be defined.** Each family needs a `#=GC RF` line; without one the run
  stops with an error unless `--rfAll` (every column with a residue) or `--rfFromOccupancy`
  (columns where at least `--rfMinOcc` of the spanning sequences have a base) is given. Either
  option replaces the RF line of every family.
- **Every AncientMode model is built with `--hand`**, including the first genome search, which
  used `--symfrac 0.3`. The starting model is built from the family's block of the prepared
  seed, so a seed file with several families also works.
- **The seed is the starting model; there is no central round.** The HMM is built from
  the seed's own match columns (`#=GC RF`, as Dfam does) and both sides extend from it,
  so every seed position stays a consensus position of the final model. This works for a
  short seed fragment and for a full model taken from another genome. The central round
  (rebuilding the model from a MAFFT alignment of the seed's copies plus flanks) kept only
  a short core for long, old repeats and could give up a side; it is removed from
  AncientMode. Both sides always extend, and a side that finds no supported extension
  stops by itself.
- **The seed Stockholm file needs a `#=GC RF` line** marking match columns with
  uppercase `X` (Dfam seeds have one). `DfamTruncationScript.sh <in.stk> <start> <end>`
  cuts a Dfam seed to a region of its consensus and writes a matching
  `.stk` and consensus `.fa`.

### Model building
- **All models are built with `hmmbuild --hand`** from match columns set by column
  occupancy (`utils/HandRF.py`, `--rfMinOcc`), instead of `--symfrac`, so HMM match
  states, consensus positions and search coordinates agree.
- **Dfam settings for every build:** `--maxinsertlen 10`, exponential entropy
  weighting (`--eentexp`, HMMER's default target of 0.62 bits/position) and HMMER's
  default sequence weights instead of `--wnone`.
- **Near-identical copies are collapsed** (segmental duplications, or one locus sampled
  twice) with cd-hit-est in the central and every recursive round
  (`utils/CollapseNearIdentical.sh`, `--maxIdentity`). With `auto`, the threshold is
  halfway between the copies' median pairwise identity and 1, within 0.80–0.98.

### Extension rounds
- **Support-based judgement of each round.** The unconserved edge of an extension
  (random flank) is trimmed by windowed conservation (`--minConservation`,
  `--supportWindow`). A round counts only if enough supported bases are added, each
  edge column backed by at least `--minCopies` copies; hits below `--minBits` are not used.
- **Retry ladder.** A failed round restarts from the last good model with a wider flank,
  more copies and a lower bit score; a side stops after `--maxRetries` failures in a
  row. Failed attempts are kept as `Failed_RoundNN_RetryN` for inspection.
- **Genome-anchored stop rule** (`utils/AnchoredSupport.py`). Each round searches the
  genome with its new bases joined to the model's outer `--anchorLength` bases. Of the
  copies spanning the `--anchorSpan` anchor bases next to the new bases, at least
  `--anchorMinRatio` of them and `--anchorMinCopies` must continue colinearly into the
  new bases (in genomes with fewer than 10 anchored copies, half of them and at least 2).
  A partly supported round keeps only its supported bases (`utils/TrimEdge.py`); the
  side stops at the second such round in a row.
- **Jump guard.** When the anchored copies rise above `--anchorMaxJump` times the
  lower count of the previous two rounds, the edge has moved into a more abundant
  neighbouring repeat: the round is discarded and the side stops.
- **Alignment cleaning is tunable.** `MSACuration.py` scores columns as
  `(1 - entropy)^a * (1 - gap fraction)^g` with a capped threshold; `a`, `g` and the cap
  are now options (`--identityExp`, `--gapExp`, `--MaxThresh`) and pipeline parameters
  for the AncientMode rounds (`--curationIdentityExp`, `--curationGapExp`,
  `--curationMaxThresh`). AncientMode now defaults to 1.4 and 0.3 (was 2 and 0.5);
  other modes keep 2 and 0.5.

### Merged model
- **Rounds are stitched, not re-aligned** (`utils/StitchMerge.py`). Re-aligning each
  round to the accumulated alignment could place the new block at the wrong end and
  duplicate part of the element. Each round's consensus is now aligned to the edge of
  the merged consensus and its match columns placed along that alignment.
- **Genome-wide refinement of the final model** (`utils/RefineModel.sh`,
  `--refineIterations`, default 8). The merged model is searched against the genome,
  the hits aligned back with hmmalign, match columns reset by occupancy and the model
  rebuilt, until the hits grow by < 2% and the consensus changes at < 3% of positions.
  An iteration that loses > 10% of the hits is discarded. Ends with few copies are
  trimmed by coverage, and blocks of a younger, abundant family the extension ran into
  are cut off (`utils/ForeignBlocks.py`). Outputs, in `Results/WorkDir/<Repeat>/Merge/`:
  `<Repeat>.Refined.hmm`, `.Consensus.fa`, `.stk` (Dfam style, match and insert
  columns), `.aln.fa`, `.RFclean.aln.fa` (consensus columns only), `.support.tsv` and
  `.iterations.tsv`.

### Examples
- `examples/AncientExtension/MamCR1.Truncated.stk` is rebuilt from the Dfam CR1_Mam seed
  (positions 1464-2204) with its `#=GC RF` line, which AncientMode needs; the previous copy
  had lost it.
- New `examples/AncientExtension/OldhAT1.Truncated.stk`: the Dfam OldhAT1 seed cut to
  positions 2000-2399, the short seed used to develop these changes (turtle genome).

### Fixes
- Extension size was overcounted when the previous consensus aligned poorly at the
  edge, so up to ~50 bases per round were stitched in twice (`utils/ExtensionSize.py`).
- The right side's search window was reset in merged-model coordinates while searching
  with the round model, so the right side always stopped at its first stall.
- A missing `-f` in the Merge step's test for the right side.
- The utils directory was resolved from inside the family directory, where it does not
  exist, so collapsing and extension measurement silently did not run.
- After a small but successful extension, the next search window is at least
  `--extension` bases wide.
