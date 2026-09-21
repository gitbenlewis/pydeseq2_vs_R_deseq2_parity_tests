# GEUVADIS R reference tables

These seven gzip-compressed TSVs preserve the R outputs from the real-Salmon
comparison introduced in parity-suite commit `b272fdc`. Decompression reproduces
the original TSV bytes; no numbers or missing values were changed.

| File | Contents |
| --- | --- |
| `r_import_counts.tsv.gz` | Unrounded gene counts from R tximport |
| `r_import_abundance.tsv.gz` | Gene abundance from R tximport |
| `r_import_length.tsv.gz` | Average transcript lengths from R tximport |
| `r_rounded_counts.tsv.gz` | DESeq2 rounded counts |
| `r_normalization_factors.tsv.gz` | Gene-by-sample normalization factors |
| `r_normalized_counts.tsv.gz` | DESeq2 normalized counts |
| `r_results.tsv.gz` | Base means and unshrunken Wald results |

The matrices have 25,343 genes and six sample columns; `r_results` has one row
per gene and the six DESeq2 result columns. `samples.tsv` defines sample order
and artificial A/B groups. This contrast is for software validation only.

`provenance.json` records source URLs and hashes, versions, analysis settings,
the tested PyDESeq2 commit, and compressed/uncompressed table checksums.
`SHA256SUMS` covers every file here except itself. No machine-specific paths,
logs, caches, or duplicate Python outputs are included.

## Use and verify

From the repository root:

```bash
(cd tests/data/geuvadis_salmon_tximport && sha256sum -c SHA256SUMS)
python -m pytest -q tests/test_geuvadis_fixtures.py
```

Pandas reads the compressed files directly:

```python
import pandas as pd
r = pd.read_csv('tests/data/geuvadis_salmon_tximport/r_results.tsv.gz',
                sep='\t', index_col='gene_id')
```

Import, rounded counts, normalization factors, normalized counts, and base means
pass the comparison gates. Downstream strict Wald gates fail: maximum LFC error
is approximately 0.0362, and adjusted-p-value NA masks disagree for 409 genes.
These reference values retain those differences. See
[`docs/salmon-parity-diagnosis.md`](../../../docs/salmon-parity-diagnosis.md).

## Regenerate deliberately

1. Use parity-suite commit `b272fdc` and the pinned `environment.yml` to reproduce
   this snapshot. Configure `config/local_env.sh` with a PyDESeq2 checkout at
   `0e56de24ee53321d52f4e85be8f410d99f764ac3`.
2. Run `bash scripts/000_run_parity.bash`. It downloads the pinned source files
   and regenerates all four comparisons. A nonzero exit is expected from the
   documented GEUVADIS Wald failures; inspect the summary and logs to distinguish
   those scientific failures from an incomplete run.
3. Preserve only the seven tables listed above from
   `results/parity/geuvadis_salmon_tximport/`. Compress each with gzip level 9
   and a zero timestamp (for example `gzip -n -9 -c INPUT.tsv > OUTPUT.tsv.gz`).
   Copy `salmon_samples.tsv` to `samples.tsv`. Do not copy `r_metadata.tsv`:
   it includes machine-specific input paths and timing data.
4. Compare decompressed table SHA-256 values with `provenance.json`. Compression
   library versions may change compressed bytes without changing table contents.
   If deliberately updating the reference snapshot, review all numerical changes
   and update versions, settings, source hashes, and table digests in provenance.
   Preserve the attribution and license; omit local paths and timings.
5. Rebuild `SHA256SUMS` for the seven archives, `samples.tsv`, `provenance.json`,
   `README.md`, and `COPYING`, then run the verification commands above. Ordinary
   parity and unit-test invocations do not update these fixtures.

## Attribution and redistribution

Source: Michael Love's tximportData, pinned at
[`d299298a82daf36f96f0506a36fa802bfddef1ff`](https://github.com/mikelove/tximportData/tree/d299298a82daf36f96f0506a36fa802bfddef1ff).
Its [DESCRIPTION](https://github.com/mikelove/tximportData/blob/d299298a82daf36f96f0506a36fa802bfddef1ff/DESCRIPTION)
declares `GPL (>= 2)`. These derived reference tables are distributed under
GPL-2.0-or-later; the GPL version 2 text is included in `COPYING`.
The original quantifications and transcript mapping remain available at the
immutable source URLs in provenance; generation code is in this repository.

The [source methods](https://github.com/mikelove/tximportData/blob/d299298a82daf36f96f0506a36fa802bfddef1ff/vignettes/tximportData.Rmd)
describe Salmon 0.6.0 quantification and the hg19 RefSeq/UCSC gene annotation.
Cite Lappalainen et al., *Transcriptome and genome sequencing uncovers functional
variation in humans*, Nature 501, 506–511 (2013), doi:10.1038/nature12531,
and the tximport and DESeq2 methods when using these outputs.

Upstream PyDESeq2 fixture integration and statistical fixes are separate tasks.
