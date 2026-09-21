"""Offline integrity and label checks for the preserved R reference outputs."""
import csv
import gzip
import hashlib
import json
from pathlib import Path


def test_geuvadis_reference_integrity_and_labels():
    directory = Path(__file__).parent / 'data' / 'geuvadis_salmon_tximport'
    checksums = {}
    for line in (directory / 'SHA256SUMS').read_text().splitlines():
        digest, name = line.split(maxsplit=1)
        assert hashlib.sha256((directory / name).read_bytes()).hexdigest() == digest
        checksums[name] = digest
    provenance = json.loads((directory / 'provenance.json').read_text())
    samples = list(csv.DictReader(
        (directory / 'samples.tsv').read_text().splitlines(), delimiter='\t',
    ))
    sample_ids = [row['sample'] for row in samples]
    assert len(sample_ids) == len(set(sample_ids)) == provenance['analysis']['expected_samples']
    assert {row['sample']: row['condition'] for row in samples} == provenance['analysis']['sample_groups']
    expected = {
        'r_import_counts', 'r_import_abundance', 'r_import_length',
        'r_rounded_counts', 'r_normalization_factors', 'r_normalized_counts', 'r_results',
    }
    assert set(provenance['tables']) == {name + '.tsv.gz' for name in expected}
    genes = None
    for name, record in provenance['tables'].items():
        raw = gzip.decompress((directory / name).read_bytes())
        assert checksums[name] == record['sha256']
        assert hashlib.sha256(raw).hexdigest() == record['uncompressed_sha256']
        rows = list(csv.reader(raw.decode().splitlines(), delimiter='\t'))
        assert rows[0][0] == 'gene_id'
        columns = ['baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj'] if name == 'r_results.tsv.gz' else sample_ids
        assert rows[0][1:] == columns
        assert all(len(row) == len(columns) + 1 for row in rows[1:])
        current = [row[0] for row in rows[1:]]
        if genes is None:
            assert len(current) == len(set(current)) == provenance['analysis']['expected_genes']
            genes = current
        assert current == genes
