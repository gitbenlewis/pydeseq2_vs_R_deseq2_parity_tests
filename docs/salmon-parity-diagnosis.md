# Real-Salmon parity diagnosis

Date: 2026-09-21. Status: diagnosis only; fixes are deferred.

## Conclusion

Keep the existing parity gates unchanged. The strongest evidence points to
differences in dispersion fitting and outlier classification, with an additional
low-count Wald standard-error difference. Salmon import and normalization pass.
The experiments below identify contributors to the failures; they do not establish
that a production fix has been implemented or that all gates would then pass.

## Scope and baseline

The comparison uses six public GEUVADIS Salmon quantifications, 25,343 genes,
and an artificial three-versus-three contrast for software validation only.
There is no biological interpretation of that contrast.

- Parity suite: commit `b272fdc`.
- PyDESeq2 checkout: `0e56de24ee53321d52f4e85be8f410d99f764ac3`.
- Python 3.11.15, R 4.5.3, DESeq2 1.50.2.
- Importers: pytximport 0.13.0 and R tximport 1.38.2.
- Source: `mikelove/tximportData@d299298a82daf36f96f0506a36fa802bfddef1ff`.
- Design: `~ condition`; B versus A; results alpha 0.1.

The canonical run passes import, rounded-count, normalization-factor,
normalized-count, and base-mean checks. Both Python input interfaces agree.
The new run fails downstream Wald gates, including maximum LFC error
approximately 0.0362 and 409 adjusted-p-value NA-mask disagreements.
The same failures persist when Python receives R-imported matrices, within
the precision of their text serialization. All three original datasets pass.

## Diagnostic evidence

### Dispersion differences explain the largest LFC discrepancy

The fitted parametric dispersion trends are approximately:

| Engine | Trend as a function of mean normalized count |
| --- | --- |
| R | `0.07570902 + 7.12299595 / mean` |
| Python | `0.079744 + 5.970422 / mean` |

The dispersion prior variances also differ: R 0.5053748 versus Python
0.4938031. Different trends and residual variances affect whether a gene is
classified as a dispersion outlier and therefore exempted from shrinkage.

| Gene | R final dispersion | Python final dispersion | R outlier | Python outlier |
| --- | ---: | ---: | --- | --- |
| GSTT2 | 2.142144 | 9.972684 | No | Yes |
| KAL1 | 2.084361 | 6.557793 | No | Yes |
| HGF | 2.671017 | 7.851784 | No | Yes |
| OR3A2 | 3.202769 | 10.000000 | No | Yes |

In a diagnostic experiment, replacing Python's final dispersions with R's
and rerunning only Python LFC fitting reduced maximum absolute LFC error
from approximately **0.0362 to 0.0000723272**. KAL1's error became
approximately 0.0000136. The maximum is below the existing 0.001 LFC gate.
This isolates final dispersion estimates as the main contributor to the
observed LFC discrepancy; substituting R estimates is not a proposed fix.

Applying R's parametric trend fitter to Python's gene-wise dispersion estimates
(finite estimates above `1e-7`) yielded coefficients approximately
`0.07639119 + 7.07827558 / mean`, much closer to R's original fit.
This implicates trend fitting in addition to differences in gene-wise estimates.
It does not isolate a single defective line or optimizer setting.

Code inspection found different iteration strategies: R initializes the trend
with `(0.1, 1)` and recomputes eligible genes from the original estimates each
iteration. Python fits before residual-based exclusion and permanently removes
excluded genes from subsequent iterations. These are candidates for controlled
follow-up experiments, not yet independently proven root causes.

Relevant PyDESeq2 methods are `DeseqDataSet._fit_parametric_dispersion_trend`
and `DeseqDataSet.fit_MAP_dispersions` in `pydeseq2/dds.py`.

### Low fitted means also affect Wald standard errors

Python's `DeseqStats.run_wald_test` computes fitted means from the coefficients
and normalization factors without applying a minimum-mean floor at that point.
For GSTT2, the minimum fitted mean is approximately 0.15739.

| GSTT2 diagnostic calculation | LFC standard error |
| --- | ---: |
| Original Python calculation | 4.205616 |
| Python with fitted means floored at 0.5 | 3.909365 |
| Python with the 0.5 floor and R's final dispersion | 2.101877 |
| R reference | 2.101893 |

The floor alone does not resolve the discrepancy. The combined substitution
closely reproduces the R standard error for this gene. OR3A2 shows similar
behavior. This supports checking minimum-mean handling alongside dispersion
fitting; it is not a full-dataset validation of a Wald fix.

### Adjusted-p-value missingness is downstream in this comparison

Applying Python's independent-filtering routine to **R's p-values and base
means** produced **zero adjusted-p-value NA-mask disagreements** with R.
The minimum retained base mean was 4.4382646509044585 in both cases.
R's reported filter cutoff was approximately 4.435323, which lies between
observed gene means.

Thus, the original 409-mask disagreement is driven by upstream input differences
to filtering in this experiment. It is not evidence, by itself, of an independent
filtering bug. This result does not claim that the filtering implementations are
equivalent on every dataset or that their adjusted values match exactly.

## Proposed follow-up — not executed

1. Compare the trend fitters using identical gene-wise dispersions and means.
   Isolate eligibility rules, initialization, convergence, and optimizer behavior.
2. Confirm whether matching the trend procedure aligns dispersion-outlier flags
   and final dispersions without importing R estimates into Python.
3. Test R-compatible minimum-mean handling in Wald standard-error calculations,
   with focused low-count examples and appropriate existing regression tests.
4. Rerun all four datasets under unchanged gates. Reassess filtering only after
   confirming the upstream p-value behavior.

Changes to PyDESeq2 require a separately authorized implementation task. Do not
relax thresholds, modify either checkout's statistical behavior, or publish fixes
as part of this diagnosis document.

## Evidence locations and limitations

The ignored directory `results/parity/geuvadis_salmon_tximport/` contains
`gate_results.tsv`, `import_comparison.json`, `shared_input_comparison.json`,
`largest_differences.tsv`, `na_mask_disagreements.tsv`, and `provenance.json`.
These are local results, not committed reference fixtures.

Diagnostic scripts and fitted objects were created in `/tmp`:
`diagnose_salmon.R`, `diagnose_salmon.py`, `salmon_counterfactual.py`,
`salmon_refit_lfc.py`, `salmon_r_stages.tsv`, `salmon_py_stages.tsv`,
`salmon_r_fit.rds`, and `salmon_py_fit.pkl`. They are temporary session artifacts
and may disappear. The measurements above preserve the findings, but this
document is not a standalone executable reproduction package.

No statistical fixes or gate changes were made during diagnosis. The complete
suite remains nonzero because the enabled real-Salmon run exposes the documented
downstream differences.
