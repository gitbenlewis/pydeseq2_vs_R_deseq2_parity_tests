#!/usr/bin/env python3
"""Run configured parity checks between PyDESeq2 and R DESeq2."""

from __future__ import annotations

import hashlib
import importlib
import importlib.metadata
import json
import logging
import os
import platform
import subprocess
import sys
import time
import urllib.request
from copy import deepcopy
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml
from scipy.stats import pearsonr
from scipy.stats import spearmanr


REPO_ROOT = Path(__file__).resolve().parent.parent
CONFIG_PATH = REPO_ROOT / "config" / "config.yaml"
with CONFIG_PATH.open(encoding="utf-8") as config_handle:
    CFG = yaml.safe_load(config_handle) or {}

PARITY_CFG = CFG["parity_params"]
RESULT_COLUMNS = (
    "baseMean",
    "log2FoldChange",
    "lfcSE",
    "stat",
    "pvalue",
    "padj",
)
THREAD_ENV_VARS = (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "BLIS_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
)


class ParityError(RuntimeError):
    """Raised when a parity input or hard gate is invalid."""


def merge_params(defaults: dict[str, Any], overrides: dict[str, Any]) -> dict[str, Any]:
    """Recursively merge one named run over shared defaults."""
    merged = deepcopy(defaults)
    for key, value in overrides.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = merge_params(merged[key], value)
        else:
            merged[key] = deepcopy(value)
    return merged


def configured_runs(parity_cfg: dict[str, Any]) -> dict[str, dict[str, Any]]:
    """Return all named runs with shared defaults applied."""
    defaults = parity_cfg.get("default_params", {})
    return {
        name: merge_params(defaults, params)
        for name, params in parity_cfg.get("task_runs", {}).items()
    }


def repo_path(value: str | Path) -> Path:
    """Resolve a config path relative to the repository."""
    path = Path(value)
    return path if path.is_absolute() else REPO_ROOT / path


def sha256_file(path: Path) -> str:
    """Return a file's SHA-256 digest."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verify_checksum(path: Path, expected: str) -> str:
    """Verify and return a file digest."""
    observed = sha256_file(path)
    if observed != expected:
        raise ParityError(
            f"SHA-256 mismatch for {path}: expected {expected}, observed {observed}"
        )
    return observed


def ensure_download(url: str, destination: Path, expected_sha256: str) -> Path:
    """Download an immutable fixture once and verify it before use."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        verify_checksum(destination, expected_sha256)
        return destination

    partial = destination.with_suffix(destination.suffix + ".part")
    try:
        urllib.request.urlretrieve(url, partial)
        verify_checksum(partial, expected_sha256)
        partial.replace(destination)
    finally:
        if partial.exists():
            partial.unlink()
    return destination


def validate_frame_alignment(
    left: pd.DataFrame,
    right: pd.DataFrame,
    label: str,
) -> None:
    """Require identical row and column identity and ordering."""
    if not left.index.equals(right.index):
        raise ParityError(f"{label} row labels or ordering differ")
    if not left.columns.equals(right.columns):
        raise ParityError(f"{label} column labels or ordering differ")


def _correlation(left: np.ndarray, right: np.ndarray, method: str) -> float:
    if np.array_equal(left, right, equal_nan=True):
        return 1.0
    if left.size < 2:
        return float("nan")
    if np.ptp(left) == 0 or np.ptp(right) == 0:
        return 1.0 if np.allclose(left, right) else float("nan")
    if method == "pearson":
        return float(pearsonr(left, right).statistic)
    return float(spearmanr(left, right).statistic)


def comparison_metrics(left: pd.DataFrame, right: pd.DataFrame) -> dict[str, Any]:
    """Summarize finite-pair differences and missing-value masks."""
    validate_frame_alignment(left, right, "comparison")
    left_values = left.to_numpy(dtype=float).ravel()
    right_values = right.to_numpy(dtype=float).ravel()
    if np.isinf(left_values).any() or np.isinf(right_values).any():
        raise ParityError("comparison contains positive or negative infinity")
    left_na = np.isnan(left_values)
    right_na = np.isnan(right_values)
    finite = np.isfinite(left_values) & np.isfinite(right_values)
    finite_left = left_values[finite]
    finite_right = right_values[finite]
    differences = np.abs(finite_right - finite_left)

    if differences.size == 0:
        quantiles = {
            "abs_diff_median": float("nan"),
            "abs_diff_p95": float("nan"),
            "abs_diff_p99": float("nan"),
            "abs_diff_max": float("nan"),
        }
    else:
        quantiles = {
            "abs_diff_median": float(np.median(differences)),
            "abs_diff_p95": float(np.quantile(differences, 0.95)),
            "abs_diff_p99": float(np.quantile(differences, 0.99)),
            "abs_diff_max": float(np.max(differences)),
        }

    return {
        "values": int(left_values.size),
        "finite_pairs": int(finite.sum()),
        "r_na": int(left_na.sum()),
        "py_na": int(right_na.sum()),
        "na_mask_disagreements": int(np.count_nonzero(left_na != right_na)),
        "pearson": _correlation(finite_left, finite_right, "pearson"),
        "spearman": _correlation(finite_left, finite_right, "spearman"),
        **quantiles,
    }


def result_decision_metrics(
    r_results: pd.DataFrame,
    py_results: pd.DataFrame,
    alphas: tuple[float, ...] = (0.05, 0.1),
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Compare LFC signs and adjusted-p-value decision sets."""
    finite = np.isfinite(r_results["log2FoldChange"]) & np.isfinite(
        py_results["log2FoldChange"]
    )
    r_sign = np.sign(r_results.loc[finite, "log2FoldChange"])
    py_sign = np.sign(py_results.loc[finite, "log2FoldChange"])
    signs = {
        "finite_pairs": int(finite.sum()),
        "concordant": int((r_sign == py_sign).sum()),
        "discordant": int((r_sign != py_sign).sum()),
        "concordance": float((r_sign == py_sign).mean()),
        "discordant_genes": list(r_results.index[finite][r_sign != py_sign]),
    }

    significant: dict[str, Any] = {}
    for alpha in alphas:
        r_set = set(
            r_results.index[
                r_results["padj"].notna() & (r_results["padj"] < alpha)
            ]
        )
        py_set = set(
            py_results.index[
                py_results["padj"].notna() & (py_results["padj"] < alpha)
            ]
        )
        union = r_set | py_set
        significant[str(alpha)] = {
            "r": len(r_set),
            "py": len(py_set),
            "overlap": len(r_set & py_set),
            "union": len(union),
            "jaccard": len(r_set & py_set) / len(union) if union else 1.0,
            "r_only": sorted(r_set - py_set),
            "py_only": sorted(py_set - r_set),
        }
    return signs, significant


def summarize_comparison(
    r_outputs: dict[str, pd.DataFrame],
    py_outputs: dict[str, pd.DataFrame],
) -> dict[str, Any]:
    """Build all metrics needed by the configured hard gates."""
    for key in ("rounded_counts", "factors", "normalized_counts", "results"):
        validate_frame_alignment(r_outputs[key], py_outputs[key], key)
    for engine, rounded in (
        ("R", r_outputs["rounded_counts"]),
        ("PyDESeq2", py_outputs["rounded_counts"]),
    ):
        if not np.isfinite(rounded.to_numpy(dtype=float)).all():
            raise ParityError(f"{engine} rounded counts contain non-finite values")

    results: dict[str, Any] = {}
    for column in RESULT_COLUMNS:
        results[column] = comparison_metrics(
            r_outputs["results"][[column]],
            py_outputs["results"][[column]],
        )
    signs, significant = result_decision_metrics(
        r_outputs["results"], py_outputs["results"]
    )
    return {
        "rounded_counts_equal": bool(
            np.array_equal(
                r_outputs["rounded_counts"].to_numpy(),
                py_outputs["rounded_counts"].to_numpy(),
                equal_nan=True,
            )
        ),
        "factors": comparison_metrics(
            r_outputs["factors"], py_outputs["factors"]
        ),
        "normalized_counts": comparison_metrics(
            r_outputs["normalized_counts"], py_outputs["normalized_counts"]
        ),
        "results": results,
        "signs": signs,
        "significant_sets": significant,
    }


def _gate(
    gates: list[dict[str, Any]],
    name: str,
    observed: Any,
    operator: str,
    threshold: Any,
) -> None:
    if operator == "min":
        passed = bool(observed >= threshold)
    elif operator == "max":
        passed = bool(observed <= threshold)
    elif operator == "equal":
        passed = bool(observed == threshold)
    else:
        raise ValueError(f"Unsupported gate operator: {operator}")
    gates.append(
        {
            "gate": name,
            "observed": observed,
            "operator": operator,
            "threshold": threshold,
            "passed": passed,
        }
    )


def _allclose_gate(
    gates: list[dict[str, Any]],
    name: str,
    left: pd.DataFrame,
    right: pd.DataFrame,
    rules: dict[str, float],
) -> None:
    validate_frame_alignment(left, right, name)
    finite = np.isfinite(left.to_numpy(dtype=float)) & np.isfinite(
        right.to_numpy(dtype=float)
    )
    passed = bool(
        finite.any()
        and np.allclose(
            left.to_numpy(dtype=float),
            right.to_numpy(dtype=float),
            rtol=float(rules["rtol"]),
            atol=float(rules["atol"]),
            equal_nan=True,
        )
    )
    max_difference = (
        float(
            np.max(
                np.abs(
                    left.to_numpy(dtype=float)[finite]
                    - right.to_numpy(dtype=float)[finite]
                )
            )
        )
        if finite.any()
        else float("nan")
    )
    gates.append(
        {
            "gate": name,
            "observed": max_difference,
            "operator": "allclose",
            "threshold": f"rtol={rules['rtol']},atol={rules['atol']}",
            "passed": passed,
        }
    )


def evaluate_gates(
    summary: dict[str, Any],
    r_outputs: dict[str, pd.DataFrame],
    py_outputs: dict[str, pd.DataFrame],
    rules: dict[str, Any],
) -> list[dict[str, Any]]:
    """Evaluate a configured gate profile without embedding thresholds in code."""
    gates: list[dict[str, Any]] = []
    _gate(
        gates,
        "rounded_counts_exact",
        summary["rounded_counts_equal"],
        "equal",
        True,
    )
    _allclose_gate(
        gates,
        "factors_allclose",
        r_outputs["factors"],
        py_outputs["factors"],
        rules["factors"],
    )
    _allclose_gate(
        gates,
        "normalized_counts_allclose",
        r_outputs["normalized_counts"],
        py_outputs["normalized_counts"],
        rules["normalized_counts"],
    )
    _allclose_gate(
        gates,
        "base_mean_allclose",
        r_outputs["results"][["baseMean"]],
        py_outputs["results"][["baseMean"]],
        rules["base_mean"],
    )

    supported_column_rules = {
        "pearson_min": ("pearson", "min"),
        "spearman_min": ("spearman", "min"),
        "abs_diff_p95_max": ("abs_diff_p95", "max"),
        "abs_diff_max": ("abs_diff_max", "max"),
        "na_mask_disagreements_max": ("na_mask_disagreements", "max"),
    }
    for column, column_rules in rules.get("result_columns", {}).items():
        for rule_name, threshold in column_rules.items():
            metric, operator = supported_column_rules[rule_name]
            _gate(
                gates,
                f"{column}_{rule_name}",
                summary["results"][column][metric],
                operator,
                threshold,
            )

    if "sign_concordance_min" in rules:
        _gate(
            gates,
            "lfc_sign_concordance",
            summary["signs"]["concordance"],
            "min",
            rules["sign_concordance_min"],
        )
    for alpha, threshold in rules.get("significant_jaccard_min", {}).items():
        _gate(
            gates,
            f"significant_jaccard_alpha_{alpha}",
            summary["significant_sets"][str(alpha)]["jaccard"],
            "min",
            threshold,
        )
    for alpha in rules.get("significant_sets_exact", []):
        decision = summary["significant_sets"][str(alpha)]
        _gate(
            gates,
            f"significant_set_exact_alpha_{alpha}",
            not decision["r_only"] and not decision["py_only"],
            "equal",
            True,
        )
    return gates


def _read_table(path: Path, index_column: str) -> pd.DataFrame:
    frame = pd.read_csv(
        path,
        sep="\t",
        index_col=index_column,
        na_values=["NA", "NaN", ""],
        keep_default_na=True,
    )
    if frame.index.has_duplicates:
        raise ParityError(f"{path} contains duplicate {index_column} values")
    return frame


def read_engine_outputs(output_dir: Path, prefix: str, mode: str) -> dict[str, pd.DataFrame]:
    """Read one engine's normalized output tables."""
    factor_name = "normalization_factors" if mode == "tximport" else "size_factors"
    factor_index = "gene_id" if mode == "tximport" else "sample"
    return {
        "rounded_counts": _read_table(
            output_dir / f"{prefix}_rounded_counts.tsv", "gene_id"
        ),
        "factors": _read_table(
            output_dir / f"{prefix}_{factor_name}.tsv", factor_index
        ),
        "normalized_counts": _read_table(
            output_dir / f"{prefix}_normalized_counts.tsv", "gene_id"
        ),
        "results": _read_table(output_dir / f"{prefix}_results.tsv", "gene_id"),
    }


def load_analysis_inputs(
    counts_path: Path,
    samples_path: Path,
    sample_column: str,
    factor: str,
    factor_levels: list[str],
    lengths_path: Path | None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame | None]:
    """Load identically ordered sample-by-gene inputs for PyDESeq2."""
    counts_table = pd.read_csv(counts_path, sep="\t")
    sample_sep = "," if samples_path.suffix == ".csv" else "\t"
    samples = pd.read_csv(samples_path, sep=sample_sep)
    if counts_table["gene_id"].duplicated().any():
        raise ParityError(f"{counts_path} contains duplicate gene identifiers")
    if samples[sample_column].duplicated().any():
        raise ParityError(f"{samples_path} contains duplicate sample identifiers")

    sample_ids = pd.Index(samples[sample_column].astype(str), name="sample")
    gene_ids = pd.Index(counts_table["gene_id"].astype(str), name="gene_id")
    missing_samples = sample_ids.difference(counts_table.columns)
    if not missing_samples.empty:
        raise ParityError(f"Counts are missing samples: {list(missing_samples)}")
    counts = counts_table.set_index("gene_id").loc[gene_ids, sample_ids].T
    counts.index = sample_ids
    counts.columns = gene_ids

    metadata = samples.set_index(sample_column).loc[sample_ids, [factor]].copy()
    metadata.index = sample_ids
    metadata[factor] = pd.Categorical(
        metadata[factor],
        categories=factor_levels,
        ordered=False,
    )
    if metadata[factor].isna().any():
        raise ParityError(
            f"Metadata factor {factor!r} contains values outside {factor_levels}"
        )

    lengths: pd.DataFrame | None = None
    if lengths_path is not None:
        lengths_table = pd.read_csv(lengths_path, sep="\t")
        length_genes = pd.Index(lengths_table["gene_id"].astype(str), name="gene_id")
        if not gene_ids.equals(length_genes):
            raise ParityError("Count and transcript-length genes are not aligned")
        missing_length_samples = sample_ids.difference(lengths_table.columns)
        if not missing_length_samples.empty:
            raise ParityError(
                f"Transcript lengths are missing samples: {list(missing_length_samples)}"
            )
        lengths = lengths_table.set_index("gene_id").loc[gene_ids, sample_ids].T
        lengths.index = sample_ids
        lengths.columns = gene_ids
    return counts, metadata, lengths


def _dense(values: Any) -> np.ndarray:
    return values.toarray() if hasattr(values, "toarray") else np.asarray(values)


def _write_frame(frame: pd.DataFrame, path: Path, index_label: str) -> None:
    frame.to_csv(
        path,
        sep="\t",
        index=True,
        index_label=index_label,
        na_rep="NA",
        float_format="%.17g",
    )


def _fit_pydeseq2(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    lengths: pd.DataFrame | None,
    params: dict[str, Any],
    *,
    use_pytximport_adata: bool,
) -> tuple[dict[str, pd.DataFrame], dict[str, Any]]:
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    start = time.perf_counter()
    if use_pytximport_adata:
        import anndata as ad

        if lengths is None:
            raise ValueError("pytximport AnnData mode requires transcript lengths")
        adata = ad.AnnData(
            X=counts.to_numpy(),
            obs=metadata.copy(),
            var=pd.DataFrame(index=counts.columns),
        )
        adata.obsm["length"] = lengths.to_numpy()
        adata.uns["counts_from_abundance"] = None
        dds = DeseqDataSet(
            adata=adata,
            design=params["design"],
            fit_type=params["fit_type"],
            size_factors_fit_type=params["size_factor_fit_type"],
            refit_cooks=params["refit_cooks"],
            n_cpus=params["n_cpus"],
            quiet=True,
        )
    else:
        dds = DeseqDataSet(
            counts=counts,
            metadata=metadata,
            transcript_lengths=lengths,
            design=params["design"],
            fit_type=params["fit_type"],
            size_factors_fit_type=params["size_factor_fit_type"],
            refit_cooks=params["refit_cooks"],
            n_cpus=params["n_cpus"],
            quiet=True,
        )
    construction_seconds = time.perf_counter() - start

    fit_start = time.perf_counter()
    dds.deseq2()
    fit_seconds = time.perf_counter() - fit_start
    observed_fit_type = str(dds.uns.get("disp_function_type"))
    validate_dispersion_fit_type(
        observed_fit_type,
        params["fit_type"],
        "PyDESeq2",
    )

    results_start = time.perf_counter()
    stats = DeseqStats(
        dds,
        contrast=[
            params["contrast_factor"],
            params["contrast_numerator"],
            params["contrast_denominator"],
        ],
        alpha=params["alpha"],
        cooks_filter=True,
        independent_filter=True,
        n_cpus=params["n_cpus"],
        quiet=True,
    )
    stats.summary()
    results_seconds = time.perf_counter() - results_start

    rounded = pd.DataFrame(
        _dense(dds.X).T,
        index=counts.columns,
        columns=counts.index,
    )
    normalized = pd.DataFrame(
        np.asarray(dds.layers["normed_counts"]).T,
        index=counts.columns,
        columns=counts.index,
    )
    if lengths is not None:
        factors = pd.DataFrame(
            np.asarray(dds.layers["normalization_factors"]).T,
            index=counts.columns,
            columns=counts.index,
        )
    else:
        factors = pd.DataFrame(
            {"size_factor": dds.obs["size_factors"].to_numpy()},
            index=counts.index,
        )
        factors.index.name = "sample"
    results = stats.results_df.loc[counts.columns, list(RESULT_COLUMNS)].copy()
    results.index.name = "gene_id"
    outputs = {
        "rounded_counts": rounded,
        "factors": factors,
        "normalized_counts": normalized,
        "results": results,
    }
    metadata_out = {
        "construction_seconds": construction_seconds,
        "fit_seconds": fit_seconds,
        "results_seconds": results_seconds,
        "total_seconds": time.perf_counter() - start,
        "genes": int(dds.n_vars),
        "samples": int(dds.n_obs),
        "design_columns": list(map(str, dds.obsm["design_matrix"].columns)),
        "dispersion_fit_type": observed_fit_type,
    }
    return outputs, metadata_out


def _assert_python_input_modes_equal(
    explicit: dict[str, pd.DataFrame],
    pytximport: dict[str, pd.DataFrame],
    rtol: float,
    atol: float,
) -> None:
    for key in explicit:
        validate_frame_alignment(explicit[key], pytximport[key], f"PyDESeq2 {key}")
        if not np.allclose(
            explicit[key].to_numpy(dtype=float),
            pytximport[key].to_numpy(dtype=float),
            rtol=rtol,
            atol=atol,
            equal_nan=True,
        ):
            raise ParityError(
                f"Explicit transcript lengths and pytximport AnnData differ for {key}"
            )


def validate_run_params(params: dict[str, Any]) -> None:
    """Reject configuration that would make the two engines scientifically differ."""
    levels = params.get("factor_levels", [])
    if len(levels) != 2 or levels[0] == levels[1]:
        raise ParityError("factor_levels must contain distinct reference and test levels")
    if params["contrast_factor"] != params["factor"]:
        raise ParityError("contrast_factor must match factor")
    if (
        params["contrast_numerator"] != levels[1]
        or params["contrast_denominator"] != levels[0]
    ):
        raise ParityError(
            "contrast must be configured test level versus reference level"
        )
    if params["fit_type"] != "parametric":
        raise ParityError("fit_type must be parametric for this parity suite")
    if params["size_factor_fit_type"] != "ratio":
        raise ParityError(
            "size_factor_fit_type must be ratio for this parity suite"
        )
    if params["n_cpus"] != 1:
        raise ParityError("n_cpus must be 1 for this parity suite")
    if not isinstance(params["refit_cooks"], bool):
        raise ParityError("refit_cooks must be a boolean")


def validate_dispersion_fit_type(
    observed: str,
    requested: str,
    engine: str,
) -> None:
    """Require the requested dispersion trend rather than an engine fallback."""
    if observed != requested:
        raise ParityError(
            f"{engine} used dispersion fit {observed!r}; requested {requested!r}"
        )


def python_runtime_versions() -> dict[str, str]:
    """Return observed Python package versions used by the comparison."""
    distributions = {
        "anndata": "anndata",
        "numpy": "numpy",
        "pandas": "pandas",
        "pydeseq2_distribution": "pydeseq2",
        "pytest": "pytest",
        "pyyaml": "PyYAML",
        "scipy": "scipy",
    }
    return {
        "python": platform.python_version(),
        **{
            name: importlib.metadata.version(distribution)
            for name, distribution in distributions.items()
        },
    }


def verify_pydeseq2_checkout(repo: Path) -> dict[str, Any]:
    """Import PyDESeq2 from the requested checkout and record Git provenance."""
    repo = repo.resolve()
    if not (repo / "pydeseq2").is_dir():
        raise ParityError(f"PYDESEQ2_REPO is not a source checkout: {repo}")
    sys.path.insert(0, str(repo))
    importlib.invalidate_caches()
    module = importlib.import_module("pydeseq2")
    module_path = Path(module.__file__).resolve()
    if not module_path.is_relative_to(repo):
        raise ParityError(
            f"Imported PyDESeq2 from {module_path}, outside requested checkout {repo}"
        )

    def git(*args: str) -> str:
        return subprocess.run(
            ["git", "-C", str(repo), *args],
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()

    return {
        "repo": str(repo),
        "module_path": str(module_path),
        "package_version": str(module.__version__),
        "git_commit": git("rev-parse", "HEAD"),
        "git_branch": git("rev-parse", "--abbrev-ref", "HEAD"),
        "git_dirty": bool(git("status", "--porcelain")),
    }


def _run_r_reference(
    params: dict[str, Any],
    paths: dict[str, Path],
    cache_dir: Path,
    output_dir: Path,
    expected_versions: dict[str, str],
    *,
    prepare_only: bool = False,
) -> None:
    script = REPO_ROOT / "scripts" / "run_r_reference.R"
    command = [
        "Rscript",
        str(script),
        "--dataset",
        params["dataset"],
        "--mode",
        params["mode"],
        "--cache-dir",
        str(cache_dir),
        "--sample-column",
        params["sample_column"],
        "--design",
        params["design"],
        "--factor",
        params["factor"],
        "--reference-level",
        params["factor_levels"][0],
        "--test-level",
        params["factor_levels"][1],
        "--contrast-factor",
        params["contrast_factor"],
        "--contrast-numerator",
        params["contrast_numerator"],
        "--contrast-denominator",
        params["contrast_denominator"],
        "--alpha",
        str(params["alpha"]),
        "--fit-type",
        params["fit_type"],
        "--size-factor-fit-type",
        params["size_factor_fit_type"],
        "--refit-cooks",
        str(params["refit_cooks"]).lower(),
        "--n-cpus",
        str(params["n_cpus"]),
        "--output-dir",
        str(output_dir),
        "--expected-r-version",
        expected_versions["r"],
        "--expected-deseq2-version",
        expected_versions["deseq2"],
        "--expected-pasilla-version",
        expected_versions["pasilla"],
        "--expected-tweedeseqcountdata-version",
        expected_versions["tweedeseqcountdata"],
        "--prepare-only",
        str(prepare_only).lower(),
    ]
    for flag, key in (
        ("--counts", "counts"),
        ("--samples", "samples"),
        ("--lengths", "lengths"),
    ):
        if key in paths:
            command.extend([flag, str(paths[key])])

    environment = os.environ.copy()
    environment.update({name: "1" for name in THREAD_ENV_VARS})
    completed = subprocess.run(
        command,
        cwd=REPO_ROOT,
        env=environment,
        capture_output=True,
        text=True,
    )
    log_prefix = "r_preparation" if prepare_only else "r"
    (output_dir / f"{log_prefix}_stdout.log").write_text(
        completed.stdout, encoding="utf-8"
    )
    (output_dir / f"{log_prefix}_stderr.log").write_text(
        completed.stderr, encoding="utf-8"
    )
    if completed.returncode:
        raise ParityError(
            f"R DESeq2 failed for {params['dataset']}: {completed.stderr[-4000:]}"
        )


def _prepare_inputs(
    params: dict[str, Any],
    dataset_cache: Path,
) -> dict[str, Path]:
    paths: dict[str, Path] = {}
    for name, spec in params.get("downloads", {}).items():
        paths[name] = ensure_download(
            spec["url"],
            dataset_cache / spec["filename"],
            spec["sha256"],
        )
    return paths


def _verify_prepared_inputs(
    params: dict[str, Any],
    dataset_cache: Path,
    paths: dict[str, Path],
) -> dict[str, Path]:
    if params["dataset"] != "srp254919":
        paths["counts"] = dataset_cache / "prepared_counts.tsv"
        paths["samples"] = dataset_cache / "prepared_samples.tsv"
    for name, expected in params.get("input_hashes", {}).items():
        verify_checksum(paths[name], expected)
    return paths


def _write_python_outputs(
    outputs: dict[str, pd.DataFrame],
    output_dir: Path,
    prefix: str,
    mode: str,
) -> None:
    _write_frame(
        outputs["rounded_counts"],
        output_dir / f"{prefix}_rounded_counts.tsv",
        "gene_id",
    )
    factor_name = "normalization_factors" if mode == "tximport" else "size_factors"
    factor_index = "gene_id" if mode == "tximport" else "sample"
    _write_frame(
        outputs["factors"],
        output_dir / f"{prefix}_{factor_name}.tsv",
        factor_index,
    )
    _write_frame(
        outputs["normalized_counts"],
        output_dir / f"{prefix}_normalized_counts.tsv",
        "gene_id",
    )
    _write_frame(
        outputs["results"],
        output_dir / f"{prefix}_results.tsv",
        "gene_id",
    )


def _largest_result_differences(
    r_results: pd.DataFrame,
    py_results: pd.DataFrame,
    count: int = 20,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for column in RESULT_COLUMNS:
        finite = np.isfinite(r_results[column]) & np.isfinite(py_results[column])
        values = pd.DataFrame(
            {
                "gene_id": r_results.index[finite],
                "r_value": r_results.loc[finite, column].to_numpy(),
                "py_value": py_results.loc[finite, column].to_numpy(),
            }
        )
        values["abs_difference"] = np.abs(values["r_value"] - values["py_value"])
        values = values.nlargest(count, "abs_difference")
        for rank, row in enumerate(values.itertuples(index=False), start=1):
            rows.append(
                {
                    "column": column,
                    "rank": rank,
                    "gene_id": row.gene_id,
                    "r_value": row.r_value,
                    "py_value": row.py_value,
                    "abs_difference": row.abs_difference,
                }
            )
    return pd.DataFrame(rows)


def _na_mask_disagreements(
    r_results: pd.DataFrame,
    py_results: pd.DataFrame,
) -> pd.DataFrame:
    """List genes whose missing-value status differs between engines."""
    rows: list[dict[str, Any]] = []
    for column in RESULT_COLUMNS:
        r_na = r_results[column].isna()
        py_na = py_results[column].isna()
        for gene_id in r_results.index[r_na != py_na]:
            rows.append(
                {
                    "column": column,
                    "gene_id": gene_id,
                    "r_is_na": bool(r_na.loc[gene_id]),
                    "py_is_na": bool(py_na.loc[gene_id]),
                    "r_value": r_results.loc[gene_id, column],
                    "py_value": py_results.loc[gene_id, column],
                }
            )
    return pd.DataFrame(
        rows,
        columns=(
            "column",
            "gene_id",
            "r_is_na",
            "py_is_na",
            "r_value",
            "py_value",
        ),
    )


def _clear_run_artifacts(output_dir: Path) -> None:
    """Remove only files owned by this runner before reusing a run directory."""
    patterns = (
        "r_*",
        "py_*",
        "comparison_summary.json",
        "gate_results.tsv",
        "known_gap.json",
        "largest_differences.tsv",
        "na_mask_disagreements.tsv",
        "provenance.json",
    )
    for pattern in patterns:
        for path in output_dir.glob(pattern):
            if path.is_file():
                path.unlink()


def _read_r_versions(output_dir: Path) -> dict[str, str]:
    metadata = pd.read_csv(output_dir / "r_metadata.tsv", sep="\t", dtype=str)
    versions = metadata.loc[
        metadata["key"].str.endswith("_version"), ["key", "value"]
    ]
    return dict(versions.itertuples(index=False, name=None))


def run_named_task(
    run_name: str,
    params: dict[str, Any],
    parity_cfg: dict[str, Any],
    pydeseq2_provenance: dict[str, Any],
) -> dict[str, Any]:
    """Run one configured R/Python comparison and persist its diagnostics."""
    output_root = repo_path(parity_cfg["output_dir"])
    cache_root = repo_path(parity_cfg["cache_dir"])
    output_dir = output_root / run_name
    dataset_cache = cache_root / params["dataset"]
    output_dir.mkdir(parents=True, exist_ok=True)
    dataset_cache.mkdir(parents=True, exist_ok=True)
    _clear_run_artifacts(output_dir)
    validate_run_params(params)

    paths = _prepare_inputs(params, dataset_cache)
    if params["dataset"] != "srp254919":
        generated_paths = (
            dataset_cache / "prepared_counts.tsv",
            dataset_cache / "prepared_samples.tsv",
        )
        if not all(path.exists() for path in generated_paths):
            _run_r_reference(
                params,
                {},
                dataset_cache,
                output_dir,
                parity_cfg["expected_versions"],
                prepare_only=True,
            )
    paths = _verify_prepared_inputs(params, dataset_cache, paths)
    _run_r_reference(
        params,
        paths,
        dataset_cache,
        output_dir,
        parity_cfg["expected_versions"],
    )

    counts, metadata, lengths = load_analysis_inputs(
        paths["counts"],
        paths["samples"],
        params["sample_column"],
        params["factor"],
        params["factor_levels"],
        paths.get("lengths"),
    )
    if counts.shape != (params["expected_samples"], params["expected_genes"]):
        raise ParityError(
            f"{run_name} dimensions are {counts.shape}, expected "
            f"({params['expected_samples']}, {params['expected_genes']})"
        )

    py_outputs, py_metadata = _fit_pydeseq2(
        counts,
        metadata,
        lengths,
        params,
        use_pytximport_adata=False,
    )
    _write_python_outputs(py_outputs, output_dir, "py", params["mode"])

    if params["mode"] == "tximport":
        adata_outputs, adata_metadata = _fit_pydeseq2(
            counts,
            metadata,
            lengths,
            params,
            use_pytximport_adata=True,
        )
        input_mode_rules = parity_cfg["pytximport_input_mode_tolerance"]
        _assert_python_input_modes_equal(
            py_outputs,
            adata_outputs,
            input_mode_rules["rtol"],
            input_mode_rules["atol"],
        )
        _write_python_outputs(
            adata_outputs, output_dir, "py_pytximport", params["mode"]
        )
        py_metadata["pytximport_adata"] = adata_metadata

    observed_python_versions = python_runtime_versions()
    py_metadata["runtime_versions"] = observed_python_versions
    py_metadata["pydeseq2_source"] = pydeseq2_provenance
    (output_dir / "py_metadata.json").write_text(
        json.dumps(py_metadata, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    r_outputs = read_engine_outputs(output_dir, "r", params["mode"])
    summary = summarize_comparison(r_outputs, py_outputs)
    gate_rules = merge_params(
        parity_cfg["gate_profiles"]["common"],
        parity_cfg["gate_profiles"][params["gate_profile"]],
    )
    gate_rules = merge_params(gate_rules, params.get("gate_overrides", {}))
    gates = evaluate_gates(summary, r_outputs, py_outputs, gate_rules)

    (output_dir / "comparison_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    pd.DataFrame(gates).to_csv(
        output_dir / "gate_results.tsv", sep="\t", index=False
    )
    _largest_result_differences(
        r_outputs["results"], py_outputs["results"]
    ).to_csv(output_dir / "largest_differences.tsv", sep="\t", index=False)
    _na_mask_disagreements(
        r_outputs["results"], py_outputs["results"]
    ).to_csv(
        output_dir / "na_mask_disagreements.tsv",
        sep="\t",
        index=False,
        na_rep="NA",
    )

    known_gap: dict[str, Any] | None = None
    if "known_gap" in params:
        gap_config = params["known_gap"]
        column = gap_config["result_column"]
        allowed = gate_rules["result_columns"][column][
            "na_mask_disagreements_max"
        ]
        known_gap = {
            **gap_config,
            "observed_na_mask_disagreements": summary["results"][column][
                "na_mask_disagreements"
            ],
            "allowed_na_mask_disagreements": allowed,
        }
        (output_dir / "known_gap.json").write_text(
            json.dumps(known_gap, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )

    input_digests = {
        name: sha256_file(path)
        for name, path in paths.items()
        if name in {"counts", "samples", "lengths"}
    }
    provenance = {
        "run_name": run_name,
        "dataset": params["dataset"],
        "config_sha256": sha256_file(CONFIG_PATH),
        "inputs": input_digests,
        "pydeseq2": pydeseq2_provenance,
        "expected_versions": parity_cfg["expected_versions"],
        "observed_python_versions": observed_python_versions,
        "observed_r_versions": _read_r_versions(output_dir),
        "known_gap": known_gap,
    }
    (output_dir / "provenance.json").write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    failed_gates = [gate["gate"] for gate in gates if not gate["passed"]]
    return {
        "run_name": run_name,
        "dataset": params["dataset"],
        "passed": not failed_gates,
        "gates": len(gates),
        "failed_gates": ",".join(failed_gates),
        "output_dir": str(output_dir),
    }


def configure_logging(output_dir: Path) -> None:
    log_dir = output_dir / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    log_path = log_dir / f"run_parity_{datetime.now():%Y%m%d_%H%M%S}.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(log_path),
            logging.StreamHandler(sys.stdout),
        ],
        force=True,
    )


def main() -> int:
    output_dir = repo_path(PARITY_CFG["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)
    summary_path = output_dir / "parity_summary.tsv"
    if summary_path.exists():
        summary_path.unlink()
    configure_logging(output_dir)
    logger = logging.getLogger(__name__)

    expected_python = PARITY_CFG["expected_versions"]["python"]
    observed_python = ".".join(map(str, sys.version_info[:3]))
    if observed_python != expected_python:
        raise ParityError(
            f"Expected Python {expected_python}, observed {observed_python}"
        )
    checkout = os.environ.get("PYDESEQ2_REPO")
    if not checkout:
        raise ParityError("PYDESEQ2_REPO must point to the checkout under test")
    pydeseq2_provenance = verify_pydeseq2_checkout(Path(checkout))
    expected_pydeseq2 = PARITY_CFG["expected_versions"]["pydeseq2"]
    if pydeseq2_provenance["package_version"] != expected_pydeseq2:
        raise ParityError(
            f"Expected PyDESeq2 {expected_pydeseq2}, observed "
            f"{pydeseq2_provenance['package_version']}"
        )
    logger.info(
        "Testing PyDESeq2 commit %s from %s",
        pydeseq2_provenance["git_commit"],
        pydeseq2_provenance["module_path"],
    )

    rows: list[dict[str, Any]] = []
    failures = 0
    for run_name, params in configured_runs(PARITY_CFG).items():
        if not params.get("run", False):
            logger.info("Skipping %s because run=false", run_name)
            rows.append(
                {
                    "run_name": run_name,
                    "dataset": params.get("dataset", ""),
                    "passed": True,
                    "gates": 0,
                    "failed_gates": "",
                    "output_dir": "",
                    "status": "skipped",
                }
            )
            continue
        logger.info("Running %s", run_name)
        try:
            row = run_named_task(
                run_name, params, PARITY_CFG, pydeseq2_provenance
            )
            row["status"] = "passed" if row["passed"] else "failed"
            failures += int(not row["passed"])
        except Exception as error:
            logger.exception("%s failed", run_name)
            failures += 1
            row = {
                "run_name": run_name,
                "dataset": params.get("dataset", ""),
                "passed": False,
                "gates": 0,
                "failed_gates": str(error),
                "output_dir": str(output_dir / run_name),
                "status": "error",
            }
        rows.append(row)

    summary = pd.DataFrame(rows)
    summary.to_csv(summary_path, sep="\t", index=False)
    if failures:
        logger.error("%d configured parity run(s) failed", failures)
        return 1
    logger.info("All configured parity runs passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
