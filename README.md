# PyDESeq2–R DESeq2 parity tests

This repository runs reproducible, config-driven parity checks between a local
PyDESeq2 source checkout and R DESeq2. It compares scientific outputs in stages:
input identity, normalization, normalized abundance, `baseMean`, and unshrunken
Wald statistics. It does not exercise the nf-core module wrappers.

The numerical-concordance checks remain the scientific acceptance gate. A
separate speed-parity runner records matched one-thread timings only after those
checks pass. Its committed calibration mode makes the reported ratios
diagnostic evidence rather than a claim that either implementation is faster.

## Datasets and provenance

All three runs are enabled by default under `parity_params.task_runs` in
`config/config.yaml`.

| Run | Pinned source and preparation | Analysis |
| --- | --- | --- |
| `srp254919_tximport` | The nf-core/test-datasets `modules` data at immutable commit `81ed58c830f2ef4640a5fd151968111dd8c5559d`. The 1,000-gene count table, six-sample sheet, and spoofed transcript-length table are downloaded to the Git-ignored cache and checked against the SHA-256 values recorded in `config/config.yaml`. | `~ treatment`; hND6 versus mCherry. R uses `DESeqDataSetFromTximport(..., countsFromAbundance = "no")`. PyDESeq2 is run both with explicit `transcript_lengths` and with `AnnData.obsm["length"]` plus `AnnData.uns["counts_from_abundance"] = None`; those Python paths must first agree with each other. |
| `pasilla` | The count matrix and sample annotation shipped in DESeq2 1.50.2. The vignette-style filter `rowSums(counts >= 10) >= 3` retains 8,148 genes across seven samples. | `~ 0 + condition`; treated versus untreated. |
| `pickrell` | `tweeDEseqCountData` 1.48.0 `pickrell.eset`. Removing genes that are zero in all samples retains 12,531 genes across 69 samples. | `~ 0 + gender`; male versus female. |

The Conda environment pins Python 3.11.15, R 4.5.3, DESeq2 1.50.2, the two
Bioconductor data packages, and the Python comparison dependencies. Generated
Pasilla and Pickrell inputs are cached outside Git. Their configured hashes are
verified before either engine analyzes the cached bytes and rechecked on every
reuse.

## Local setup and canonical run

Create the pinned comparison environment:

```bash
conda env create --file environment.yml
```

Create the machine-local bootstrap file and edit its three values:

```bash
cp config/local_env.sh.example config/local_env.sh
```

`config/local_env.sh` is ignored by Git. Set:

```bash
export PYDESEQ2_REPO=/home/ubuntu/PyDESeq2-pytximport
export PARITY_CONDA_ENV=nfcore-pydeseq2-comparison
export CONDA_BASE=/home/ubuntu/miniconda3
```

`PYDESEQ2_REPO` may point at any checkout to test, but it must contain the
PyDESeq2 source tree and Git metadata. From this repository's root, run:

```bash
bash scripts/000_run_parity.bash
```

This is the canonical entrypoint used by CI. It activates
`PARITY_CONDA_ENV`, limits native math libraries to one thread, places
`PYDESEQ2_REPO` first on `PYTHONPATH`, runs this repository's focused unit
tests and the checkout's transcript-length normalization tests, and then
executes all enabled parity runs. The Python runner resolves
`pydeseq2.__file__` and fails unless it is inside the resolved
`PYDESEQ2_REPO`; each report also records the checkout's actual Git SHA. This
prevents an installed PyDESeq2 package from silently replacing the requested
source checkout.

The Git-ignored cache defaults to `.cache/parity` and results default to
`results/parity`. Both locations are configurable in `config/config.yaml`.
The first run needs network access for the pinned SRP254919 files; later runs
reuse the verified cache.

## Speed-parity benchmark

Run the canonical speed suite with:

```bash
bash scripts/010_run_speed_parity.bash
```

This entrypoint first runs `scripts/000_run_parity.bash`, so the benchmark is
accepted only when fresh numerical-parity outputs match the current config,
input hashes, and PyDESeq2 source SHA. It then benchmarks all three datasets
with one unmeasured warm-up and seven measured repetitions. Every R or Python
cell runs in a fresh process, execution order is deterministically
counterbalanced, analyses are serialized, and native math libraries remain
limited to one thread.

The primary timing boundary is `fit_seconds + results_seconds`: R times
`DESeq()` plus `results()`, while Python times `DeseqDataSet.deseq2()` plus
`DeseqStats.summary()`. Construction time and fresh-process wall time are
recorded separately; input preparation, interpreter startup, and diagnostic
file writing are not part of the primary comparison. The paired ratio is
`exp(median(log(PyDESeq2 core / R core)))` across aligned repetitions.
Fresh-process wall time remains an engine-local startup diagnostic and is not
used to compare implementations or enforce a ratio.

For `srp254919_tximport`, explicit `transcript_lengths` is the primary Python
cell compared with R's tximport path. The compatible AnnData input path is a
separate diagnostic cell, reported as AnnData/explicit Python overhead; it is
not combined with the primary ratio. Pasilla and Pickrell each compare their
matrix-based Python analysis directly with R.

The committed `calibration_mode: true` reports, but does not enforce, the
provisional PyDESeq2/R guides of 1.50 for SRP254919, 6.00 for Pasilla, and 2.75
for Pickrell, the provisional SRP AnnData/explicit guide of 1.15, and the
aspirational cross-dataset ratio of 1.25. These values are not calibrated
acceptance thresholds. Hard performance limits require repeated baseline
sessions on a stable runner: collect at least 20 clean invocations across five
days before freezing dataset-specific limits from the robust distribution.
Hard limits are intentionally disabled on shared GitHub hosts.

No timing observation is discarded. The runner reports medians, ranges, IQRs,
MADs, and robust coefficients of variation for construction, fit, results,
core, and process-wall timings. If any core-timing robust CV exceeds 0.15, it
repeats the complete dataset once. A second noisy attempt is classified as
`inconclusive_infrastructure`, not as an implementation regression.

## Continuous integration

Correctness CI checks out `gitbenlewis/PyDESeq2@main` beside this repository,
verifies that all three configured runs remain enabled, and invokes
`scripts/000_run_parity.bash`. That entrypoint includes the checkout's focused
`test_transcript_length_normalization.py` suite (current baseline: 41 passed
and 12 skipped).

The speed workflow uses the same source checkout, pinned Conda environment,
verified input cache, and `scripts/010_run_speed_parity.bash` entrypoint. Pull
requests use a shorter one-warm-up/three-measurement diagnostic invocation;
scheduled and manual runs use the YAML-backed canonical one-plus-seven
configuration. Because GitHub-hosted runners are variable, calibration mode
remains required and CI makes no hard ratio or speed-superiority claim.
Correctness and speed diagnostics are uploaded whether the run passes or
fails.

## Configuration and run flags

`config/config.yaml` contains one `parity_params` block:

- `output_dir` and `cache_dir` select the untracked result and input-cache
  roots; `expected_versions` pins the analysis engines and data packages.
- `default_params` defines shared fitting behavior,
  `pytximport_input_mode_tolerance` controls agreement between the two Python
  input paths, and `gate_profiles` defines the R-versus-Python thresholds.
- `task_runs` contains `srp254919_tximport`, `pasilla`, and `pickrell`.
- Each named run merges its values over `default_params`. Set its `run` value
  to `false` only for a focused local iteration; the committed defaults and CI
  keep all three values `true`.
- `speed` defines the benchmark repetitions, counterbalancing seed, one-thread
  contract, calibration guides, noise policy, output root, and per-dataset run
  flags. Its named runs reference the corresponding scientific runs rather
  than duplicating dataset, design, contrast, or hash settings.

Dataset paths, expected dimensions, source hashes, designs, contrasts, and
acceptance thresholds belong in the YAML rather than in the Bash entrypoint.
The runner logs a skipped run when `run: false`.

## Outputs and hard gates

`results/parity/parity_summary.tsv` gives the status of every configured run.
Each `results/parity/<run>/` directory contains:

- `r_*` and `py_*` TSVs for rounded counts, normalization factors or size
  factors, normalized counts, and results. SRP254919 also has
  `py_pytximport_*` TSVs for the compatible AnnData input path.
- `comparison_summary.json`, `gate_results.tsv`, and
  `largest_differences.tsv`; `na_mask_disagreements.tsv` names every gene with
  a missing-value-mask disagreement.
- `provenance.json`, `py_metadata.json`, `r_metadata.tsv`, and captured R
  session information, standard-output, and error logs. Pickrell also writes
  `known_gap.json` with the observed and permitted adjusted-p-value mask gap.

Together these diagnostics record input and config hashes, expected and
observed engine versions, the PyDESeq2 source path and Git SHA, timing metadata,
and the largest numerical disagreements. They are intentionally untracked.
Aggregate runner logs are written below `results/parity/logs/`, with the
canonical Bash transcript below `scripts/logs/`.

Each speed invocation is immutable below
`results/speed_parity/<UTC invocation ID>/`; `latest_run.json` points to the
newest bundle. `speed_trials.tsv` contains every warm-up and measured cell with
its execution order, phase timings, process wall time, dimensions, scientific
output fingerprint, and worker directory. `speed_attempts.tsv` contains the
robust summaries and ratio observations for every attempt, while
`speed_summary.tsv` and `speed_summary.json` record the selected attempt and
final status for all three datasets. Per-worker logs and metadata remain under
the dataset attempt directories.

The speed `provenance.json` records the config hash, command, PyDESeq2 source
path and Git SHA, Python package versions, host platform and CPU, CPU affinity,
BLAS details, native thread settings, and the exact passing correctness
provenance used as its precondition. The invocation log is stored alongside
these files, and the Bash transcript remains below `scripts/logs/`. All speed
artifacts are intentionally untracked.

The suite fails immediately on different sample or gene labels, dimensions, or
rounded counts, or if either engine falls back from the configured parametric
dispersion trend. Normalization factors (or scalar size factors) use
`rtol = atol = 1e-12`; normalized counts and `baseMean` use `rtol = 1e-10` and
`atol = 1e-8`.

For `srp254919_tximport`, p-value and adjusted-p-value NA masks, LFC signs, and
significant sets at alpha 0.05 and 0.1 must match exactly. Its configured gates
also require:

| Quantity | Correlation gate | Error gate |
| --- | --- | --- |
| LFC | Pearson ≥ 0.999999 | max ≤ 0.001 |
| LFC SE | Pearson ≥ 0.9999 | max ≤ 0.02 |
| Wald statistic | Pearson ≥ 0.9999 | max ≤ 0.05 |
| p-value | Spearman ≥ 0.9999 | max ≤ 0.005 |
| adjusted p-value | Spearman ≥ 0.9999 | max ≤ 0.01 |

For Pasilla and Pickrell, LFC Pearson and Spearman correlations must be at least
0.9999, with p95/max absolute errors no greater than 0.0005/0.025 and no LFC NA
mask disagreements. P-value Pearson and Spearman correlations must be at least
0.9997, with p95/max errors no greater than 0.02/0.20 and at most one NA-mask
disagreement. Adjusted-p-value Pearson and Spearman correlations must be at
least 0.98/0.998, with p95/max errors no greater than 0.10/0.25. LFC sign
concordance must be at least 0.999 and the significant-set Jaccard index at
alpha 0.1 must be at least 0.97.

Pasilla permits no adjusted-p-value NA-mask disagreements. Pickrell explicitly
permits at most 1,250 as a guard for the known independent-filtering difference
that can affect part of this dataset. `known_gap.json` records the observed
count beside that allowance on every run. This does not relax the other
Pickrell gates or imply exact adjusted-p-value parity.

## Scope

The speed suite intentionally excludes calibrated hard performance assertions,
memory/RSS gates, multi-CPU scaling, and Nextflow/module-wrapper overhead. The
scientific suite still excludes VST, coefficient shrinkage, likelihood-ratio
tests, and blocked or control-gene normalization variants. Neither suite
changes the nf-core/modules checkout or the PyDESeq2 checkout under test.
