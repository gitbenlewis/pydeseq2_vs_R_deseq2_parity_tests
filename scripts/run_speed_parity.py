#!/usr/bin/env python3
"""Benchmark equivalent PyDESeq2 and R DESeq2 analysis phases."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import math
import os
import platform
import subprocess
import sys
import time
from copy import deepcopy
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable

import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))
import run_parity  # noqa: E402


REPO_ROOT = SCRIPT_DIR.parent
PARITY_CFG = run_parity.PARITY_CFG
SPEED_CFG = PARITY_CFG["speed"]
THREAD_ENV_VARS = run_parity.THREAD_ENV_VARS
TIMING_FIELDS = (
    "construction_seconds",
    "fit_seconds",
    "results_seconds",
)
RAW_COLUMNS = (
    "run_name",
    "dataset",
    "engine",
    "python_input_mode",
    "attempt",
    "phase",
    "repetition",
    "order_position",
    "execution_sequence",
    "construction_seconds",
    "fit_seconds",
    "results_seconds",
    "core_seconds",
    "process_wall_seconds",
    "genes",
    "samples",
    "dispersion_fit_type",
    "fingerprint",
    "worker_output_dir",
)


class SpeedParityError(run_parity.ParityError):
    """Raised when benchmark evidence is invalid or incomplete."""


def _positive_float(value: Any, label: str) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError) as error:
        raise SpeedParityError(f"{label} must be numeric") from error
    if not math.isfinite(parsed) or parsed <= 0:
        raise SpeedParityError(f"{label} must be finite and greater than zero")
    return parsed


def _nonnegative_int(value: Any, label: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise SpeedParityError(f"{label} must be a non-negative integer")
    return value


def validate_speed_config(speed_cfg: dict[str, Any]) -> None:
    """Validate execution settings that affect benchmark comparability."""
    if not isinstance(speed_cfg.get("run"), bool):
        raise SpeedParityError("speed.run must be a boolean")
    _nonnegative_int(speed_cfg.get("warmup_repetitions"), "warmup_repetitions")
    measured = _nonnegative_int(
        speed_cfg.get("measured_repetitions"), "measured_repetitions"
    )
    if measured == 0:
        raise SpeedParityError("measured_repetitions must be greater than zero")
    if speed_cfg.get("execution_order") != "counterbalanced":
        raise SpeedParityError("execution_order must be counterbalanced")
    _nonnegative_int(speed_cfg.get("random_seed"), "random_seed")
    if speed_cfg.get("n_cpus") != 1:
        raise SpeedParityError("speed.n_cpus must be 1")
    if speed_cfg.get("native_math_threads") != 1:
        raise SpeedParityError("native_math_threads must be 1")
    if not isinstance(speed_cfg.get("calibration_mode"), bool):
        raise SpeedParityError("calibration_mode must be a boolean")
    _positive_float(
        speed_cfg.get("aspirational_ratio_max"), "aspirational_ratio_max"
    )
    noise = speed_cfg.get("noise", {})
    max_robust_cv = _positive_float(
        noise.get("max_robust_cv"), "noise.max_robust_cv"
    )
    if max_robust_cv >= 1:
        raise SpeedParityError("noise.max_robust_cv must be less than one")
    if not isinstance(noise.get("retry_noisy_once"), bool):
        raise SpeedParityError("noise.retry_noisy_once must be a boolean")
    if not isinstance(speed_cfg.get("default_params", {}).get("run"), bool):
        raise SpeedParityError("speed.default_params.run must be a boolean")
    for run_name, params in speed_cfg.get("task_runs", {}).items():
        if not isinstance(params.get("run", True), bool):
            raise SpeedParityError(f"{run_name}.run must be a boolean")
        if not isinstance(params.get("benchmark_adata", False), bool):
            raise SpeedParityError(f"{run_name}.benchmark_adata must be a boolean")
        _positive_float(
            params.get("provisional_ratio_max"),
            f"{run_name}.provisional_ratio_max",
        )
        if params.get("benchmark_adata", False):
            _positive_float(
                params.get("provisional_adata_overhead_max"),
                f"{run_name}.provisional_adata_overhead_max",
            )


def configured_speed_runs(
    parity_cfg: dict[str, Any],
) -> dict[str, dict[str, dict[str, Any]]]:
    """Return speed settings paired with the existing scientific run settings."""
    speed_cfg = parity_cfg["speed"]
    validate_speed_config(speed_cfg)
    scientific_runs = run_parity.configured_runs(parity_cfg)
    speed_defaults = speed_cfg.get("default_params", {})
    configured: dict[str, dict[str, dict[str, Any]]] = {}
    for run_name, speed_overrides in speed_cfg.get("task_runs", {}).items():
        if run_name not in scientific_runs:
            raise SpeedParityError(
                f"Speed run {run_name!r} has no matching parity task"
            )
        speed_params = run_parity.merge_params(speed_defaults, speed_overrides)
        speed_params["run"] = bool(
            speed_cfg["run"]
            and scientific_runs[run_name].get("run", False)
            and speed_params.get("run", False)
        )
        configured[run_name] = {
            "analysis": scientific_runs[run_name],
            "speed": speed_params,
        }
    return configured


def validate_thread_environment(expected: int = 1) -> dict[str, str]:
    """Require every configured native thread cap to match the benchmark."""
    observed = {name: os.environ.get(name, "") for name in THREAD_ENV_VARS}
    mismatches = {
        name: value for name, value in observed.items() if value != str(expected)
    }
    if mismatches:
        raise SpeedParityError(
            f"Native thread settings must all equal {expected}: {mismatches}"
        )
    return observed


def validate_python_environment(
    source: dict[str, Any],
    parity_cfg: dict[str, Any],
) -> dict[str, str]:
    """Require the same pinned Python and PyDESeq2 versions as correctness."""
    versions = run_parity.python_runtime_versions()
    expected = parity_cfg["expected_versions"]
    if versions["python"] != expected["python"]:
        raise SpeedParityError(
            f"Expected Python {expected['python']}, observed {versions['python']}"
        )
    if source["package_version"] != expected["pydeseq2"]:
        raise SpeedParityError(
            f"Expected PyDESeq2 {expected['pydeseq2']}, observed "
            f"{source['package_version']}"
        )
    return versions


def _stable_offset(run_name: str, random_seed: int) -> int:
    digest = hashlib.sha256(run_name.encode("utf-8")).hexdigest()
    return int(digest[:8], 16) + int(random_seed)


def build_schedule(
    run_name: str,
    warmup_repetitions: int,
    measured_repetitions: int,
    attempt: int,
    benchmark_adata: bool,
    random_seed: int,
) -> list[dict[str, Any]]:
    """Build a deterministic rotation of fresh-process implementation cells."""
    cells = [
        ("r_deseq2", ""),
        ("pydeseq2", "explicit"),
    ]
    if benchmark_adata:
        cells.append(("pydeseq2", "anndata"))
    schedule: list[dict[str, Any]] = []
    sequence = 0
    base_offset = _stable_offset(run_name, random_seed) + attempt - 1
    for phase, repetitions in (
        ("warmup", warmup_repetitions),
        ("measured", measured_repetitions),
    ):
        phase_offset = 0 if phase == "warmup" else warmup_repetitions
        for repetition in range(1, repetitions + 1):
            offset = (base_offset + phase_offset + repetition - 1) % len(cells)
            ordered = cells[offset:] + cells[:offset]
            for order_position, (engine, input_mode) in enumerate(ordered, start=1):
                sequence += 1
                schedule.append(
                    {
                        "engine": engine,
                        "python_input_mode": input_mode,
                        "attempt": attempt,
                        "phase": phase,
                        "repetition": repetition,
                        "order_position": order_position,
                        "execution_sequence": sequence,
                    }
                )
    return schedule


def _analysis_paths(
    params: dict[str, Any],
) -> tuple[dict[str, Path], Path]:
    dataset_cache = (
        run_parity.repo_path(PARITY_CFG["cache_dir"]) / params["dataset"]
    )
    paths = run_parity._prepare_inputs(params, dataset_cache)
    paths = run_parity._verify_prepared_inputs(params, dataset_cache, paths)
    return paths, dataset_cache


def _scientific_output_paths(
    output_dir: Path,
    prefix: str,
    mode: str,
) -> list[Path]:
    factor_name = "normalization_factors" if mode == "tximport" else "size_factors"
    return [
        output_dir / f"{prefix}_rounded_counts.tsv",
        output_dir / f"{prefix}_{factor_name}.tsv",
        output_dir / f"{prefix}_normalized_counts.tsv",
        output_dir / f"{prefix}_results.tsv",
    ]


def fingerprint_outputs(paths: list[Path], *, remove: bool = False) -> str:
    """Hash only deterministic scientific outputs, optionally removing them."""
    digest = hashlib.sha256()
    for position, path in enumerate(paths):
        if not path.is_file():
            raise SpeedParityError(f"Missing scientific output: {path}")
        digest.update(str(position).encode("ascii"))
        digest.update(run_parity.sha256_file(path).encode("ascii"))
    fingerprint = digest.hexdigest()
    if remove:
        for path in paths:
            path.unlink()
    return fingerprint


def parse_r_metadata(
    path: Path,
    params: dict[str, Any],
) -> dict[str, Any]:
    """Read and validate one R worker's timing and scientific identity."""
    frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    if list(frame.columns) != ["key", "value"]:
        raise SpeedParityError(f"Unexpected R metadata columns in {path}")
    if frame["key"].duplicated().any():
        raise SpeedParityError(f"Duplicate R metadata keys in {path}")
    metadata = dict(frame.itertuples(index=False, name=None))
    expected = {
        "dataset": str(params["dataset"]),
        "mode": str(params["mode"]),
        "genes": str(params["expected_genes"]),
        "samples": str(params["expected_samples"]),
        "design": str(params["design"]),
        "contrast": (
            f"{params['contrast_factor']}_{params['contrast_numerator']}"
            f"_vs_{params['contrast_denominator']}"
        ),
        "dispersion_fit_type": str(params["fit_type"]),
        "size_factor_type": str(params["size_factor_fit_type"]),
        "n_cpus": "1",
        "parallel": "FALSE",
    }
    for key, expected_value in expected.items():
        if metadata.get(key) != expected_value:
            raise SpeedParityError(
                f"R metadata {key}={metadata.get(key)!r}, "
                f"expected {expected_value!r}"
            )
    timings = {
        key: _positive_float(metadata.get(key), f"R {key}")
        for key in TIMING_FIELDS
    }
    return {
        **timings,
        "genes": int(metadata["genes"]),
        "samples": int(metadata["samples"]),
        "dispersion_fit_type": metadata["dispersion_fit_type"],
    }


def timing_record(
    metadata: dict[str, Any],
    *,
    run_name: str,
    dataset: str,
    engine: str,
    input_mode: str,
    process_wall_seconds: float,
    fingerprint: str,
    worker_output_dir: Path,
) -> dict[str, Any]:
    """Normalize one engine's phase timings into the raw benchmark schema."""
    timings = {
        key: _positive_float(metadata.get(key), f"{engine} {key}")
        for key in TIMING_FIELDS
    }
    process_wall = _positive_float(process_wall_seconds, "process_wall_seconds")
    genes = int(metadata["genes"])
    samples = int(metadata["samples"])
    return {
        "run_name": run_name,
        "dataset": dataset,
        "engine": engine,
        "python_input_mode": input_mode,
        **timings,
        "core_seconds": timings["fit_seconds"] + timings["results_seconds"],
        "process_wall_seconds": process_wall,
        "genes": genes,
        "samples": samples,
        "dispersion_fit_type": str(metadata["dispersion_fit_type"]),
        "fingerprint": fingerprint,
        "worker_output_dir": str(worker_output_dir),
    }


def run_python_worker(run_name: str, input_mode: str, output_dir: Path) -> int:
    """Execute one fresh-process PyDESeq2 analysis and persist compact metadata."""
    if input_mode not in {"explicit", "anndata"}:
        raise SpeedParityError(f"Unsupported Python input mode: {input_mode}")
    runs = configured_speed_runs(PARITY_CFG)
    if run_name not in runs:
        raise SpeedParityError(f"Unknown speed run: {run_name}")
    params = runs[run_name]["analysis"]
    run_parity.validate_run_params(params)
    validate_thread_environment(SPEED_CFG["native_math_threads"])
    checkout = os.environ.get("PYDESEQ2_REPO")
    if not checkout:
        raise SpeedParityError("PYDESEQ2_REPO must point to the checkout under test")
    source = run_parity.verify_pydeseq2_checkout(Path(checkout))
    validate_python_environment(source, PARITY_CFG)
    paths, _ = _analysis_paths(params)
    counts, metadata, lengths = run_parity.load_analysis_inputs(
        paths["counts"],
        paths["samples"],
        params["sample_column"],
        params["factor"],
        params["factor_levels"],
        paths.get("lengths"),
    )
    expected_shape = (params["expected_samples"], params["expected_genes"])
    if counts.shape != expected_shape:
        raise SpeedParityError(
            f"{run_name} dimensions are {counts.shape}, expected {expected_shape}"
        )
    if input_mode == "anndata" and params["mode"] != "tximport":
        raise SpeedParityError("AnnData timing is only valid for tximport runs")

    output_dir.mkdir(parents=True, exist_ok=False)
    outputs, timing = run_parity._fit_pydeseq2(
        counts,
        metadata,
        lengths,
        params,
        use_pytximport_adata=input_mode == "anndata",
    )
    prefix = "py"
    run_parity._write_python_outputs(outputs, output_dir, prefix, params["mode"])
    fingerprint = fingerprint_outputs(
        _scientific_output_paths(output_dir, prefix, params["mode"]),
        remove=True,
    )
    worker_metadata = {
        **timing,
        "fingerprint": fingerprint,
        "input_mode": input_mode,
        "pydeseq2_source": source,
    }
    (output_dir / "worker_metadata.json").write_text(
        json.dumps(worker_metadata, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return 0


def _run_python_trial(
    run_name: str,
    params: dict[str, Any],
    input_mode: str,
    output_dir: Path,
    expected_source: dict[str, Any],
) -> dict[str, Any]:
    command = [
        sys.executable,
        str(Path(__file__).resolve()),
        "--worker",
        "--run-name",
        run_name,
        "--input-mode",
        input_mode,
        "--worker-output-dir",
        str(output_dir),
    ]
    started = time.perf_counter()
    completed = subprocess.run(
        command,
        cwd=REPO_ROOT,
        env=os.environ.copy(),
        capture_output=True,
        text=True,
    )
    process_wall = time.perf_counter() - started
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    (output_dir.parent / f"{output_dir.name}_stdout.log").write_text(
        completed.stdout, encoding="utf-8"
    )
    (output_dir.parent / f"{output_dir.name}_stderr.log").write_text(
        completed.stderr, encoding="utf-8"
    )
    if completed.returncode:
        raise SpeedParityError(
            f"PyDESeq2 worker failed for {run_name}/{input_mode}: "
            f"{completed.stderr[-4000:]}"
        )
    metadata = json.loads(
        (output_dir / "worker_metadata.json").read_text(encoding="utf-8")
    )
    if metadata.get("input_mode") != input_mode:
        raise SpeedParityError("Python worker input mode changed")
    if (
        int(metadata.get("genes", -1)) != params["expected_genes"]
        or int(metadata.get("samples", -1)) != params["expected_samples"]
    ):
        raise SpeedParityError("Python worker dimensions changed")
    if metadata.get("dispersion_fit_type") != params["fit_type"]:
        raise SpeedParityError("Python worker dispersion fit type changed")
    worker_source = metadata.get("pydeseq2_source", {})
    source_keys = ("git_commit", "module_path", "package_version", "git_dirty")
    if any(
        worker_source.get(key) != expected_source.get(key)
        for key in source_keys
    ):
        raise SpeedParityError("PyDESeq2 source changed during the benchmark")
    return timing_record(
        metadata,
        run_name=run_name,
        dataset=params["dataset"],
        engine="pydeseq2",
        input_mode=input_mode,
        process_wall_seconds=process_wall,
        fingerprint=str(metadata["fingerprint"]),
        worker_output_dir=output_dir,
    )


def _run_r_trial(
    run_name: str,
    params: dict[str, Any],
    output_dir: Path,
) -> dict[str, Any]:
    started = time.perf_counter()
    paths, dataset_cache = _analysis_paths(params)
    output_dir.mkdir(parents=True, exist_ok=False)
    run_parity._run_r_reference(
        params,
        paths,
        dataset_cache,
        output_dir,
        PARITY_CFG["expected_versions"],
    )
    metadata = parse_r_metadata(output_dir / "r_metadata.tsv", params)
    fingerprint = fingerprint_outputs(
        _scientific_output_paths(output_dir, "r", params["mode"]),
        remove=True,
    )
    process_wall = time.perf_counter() - started
    return timing_record(
        metadata,
        run_name=run_name,
        dataset=params["dataset"],
        engine="r_deseq2",
        input_mode="",
        process_wall_seconds=process_wall,
        fingerprint=fingerprint,
        worker_output_dir=output_dir,
    )


def run_attempt(
    run_name: str,
    params: dict[str, Any],
    speed_params: dict[str, Any],
    speed_cfg: dict[str, Any],
    attempt: int,
    output_dir: Path,
    trial_runner: Callable[..., dict[str, Any]] | None = None,
    record_callback: Callable[[dict[str, Any]], None] | None = None,
    expected_source: dict[str, Any] | None = None,
) -> list[dict[str, Any]]:
    """Execute one complete warm-up and measured attempt."""
    schedule = build_schedule(
        run_name,
        speed_cfg["warmup_repetitions"],
        speed_cfg["measured_repetitions"],
        attempt,
        bool(speed_params.get("benchmark_adata", False)),
        speed_cfg["random_seed"],
    )
    records: list[dict[str, Any]] = []
    for item in schedule:
        label = (
            f"{item['phase']}_{item['repetition']:03d}_"
            f"{item['order_position']:02d}_{item['engine']}"
        )
        if item["python_input_mode"]:
            label += f"_{item['python_input_mode']}"
        worker_dir = output_dir / f"attempt_{attempt}" / label
        if trial_runner is not None:
            record = trial_runner(
                run_name,
                params,
                item["engine"],
                item["python_input_mode"],
                worker_dir,
            )
        elif item["engine"] == "r_deseq2":
            record = _run_r_trial(run_name, params, worker_dir)
        else:
            if expected_source is None:
                raise SpeedParityError(
                    "Real Python trials require expected source provenance"
                )
            record = _run_python_trial(
                run_name,
                params,
                item["python_input_mode"],
                worker_dir,
                expected_source,
            )
        completed_record = {**record, **item}
        records.append(completed_record)
        if record_callback is not None:
            record_callback(completed_record)
    return records


def validate_trial_records(
    records: list[dict[str, Any]],
    run_name: str,
    params: dict[str, Any],
    speed_params: dict[str, Any],
    speed_cfg: dict[str, Any],
    attempt: int,
    expected_fingerprints: dict[str, str] | None = None,
) -> pd.DataFrame:
    """Reject incomplete, duplicated, nondeterministic, or invalid trials."""
    frame = pd.DataFrame(records)
    missing_columns = set(RAW_COLUMNS).difference(frame.columns)
    if missing_columns:
        raise SpeedParityError(
            f"{run_name} trial records lack columns: {sorted(missing_columns)}"
        )
    key_columns = [
        "attempt",
        "phase",
        "repetition",
        "engine",
        "python_input_mode",
    ]
    if frame.duplicated(key_columns).any():
        raise SpeedParityError(f"{run_name} contains duplicate benchmark cells")
    expected = build_schedule(
        run_name,
        speed_cfg["warmup_repetitions"],
        speed_cfg["measured_repetitions"],
        attempt,
        bool(speed_params.get("benchmark_adata", False)),
        speed_cfg["random_seed"],
    )
    expected_cells = {
        (
            item["attempt"],
            item["phase"],
            item["repetition"],
            item["engine"],
            item["python_input_mode"],
        )
        for item in expected
    }
    observed_cells = {
        tuple(row)
        for row in frame[key_columns].itertuples(index=False, name=None)
    }
    if observed_cells != expected_cells:
        raise SpeedParityError(f"{run_name} has incomplete benchmark cells")
    for column in (*TIMING_FIELDS, "core_seconds", "process_wall_seconds"):
        values = pd.to_numeric(frame[column], errors="coerce").to_numpy()
        if not np.isfinite(values).all() or np.any(values <= 0):
            raise SpeedParityError(f"{run_name} contains invalid {column}")
    if not (frame["genes"].astype(int) == params["expected_genes"]).all():
        raise SpeedParityError(f"{run_name} gene dimensions changed")
    if not (frame["samples"].astype(int) == params["expected_samples"]).all():
        raise SpeedParityError(f"{run_name} sample dimensions changed")
    if not (frame["dispersion_fit_type"] == params["fit_type"]).all():
        raise SpeedParityError(f"{run_name} dispersion fit type changed")
    for (engine, input_mode), group in frame.groupby(
        ["engine", "python_input_mode"], dropna=False
    ):
        if group["fingerprint"].nunique(dropna=False) != 1:
            raise SpeedParityError(
                f"{run_name} produced nondeterministic scientific outputs"
            )
        if expected_fingerprints is not None:
            cell = (
                "r_deseq2"
                if engine == "r_deseq2"
                else f"pydeseq2_{input_mode}"
            )
            if cell not in expected_fingerprints:
                raise SpeedParityError(
                    f"{run_name} lacks a correctness fingerprint for {cell}"
                )
            if group["fingerprint"].iloc[0] != expected_fingerprints[cell]:
                raise SpeedParityError(
                    f"{run_name} {cell} differs from numerical parity output"
                )
    return frame


def robust_timing_summary(values: np.ndarray | list[float]) -> dict[str, float]:
    """Summarize positive durations without dropping observations."""
    array = np.asarray(values, dtype=float)
    if (
        array.size == 0
        or not np.isfinite(array).all()
        or np.any(array <= 0)
    ):
        raise SpeedParityError("Timing summary requires positive finite values")
    median = float(np.median(array))
    mad = float(np.median(np.abs(array - median)))
    q1, q3 = np.quantile(array, [0.25, 0.75])
    return {
        "median": median,
        "minimum": float(np.min(array)),
        "maximum": float(np.max(array)),
        "iqr": float(q3 - q1),
        "mad": mad,
        "robust_cv": float(1.4826 * mad / median),
    }


def paired_median_ratio(
    numerator: pd.Series,
    denominator: pd.Series,
) -> float:
    """Return exp(median(log(numerator / denominator))) for paired trials."""
    if not numerator.index.equals(denominator.index):
        raise SpeedParityError("Paired timing repetitions are not aligned")
    numerator_values = numerator.to_numpy(dtype=float)
    denominator_values = denominator.to_numpy(dtype=float)
    if (
        not np.isfinite(numerator_values).all()
        or not np.isfinite(denominator_values).all()
        or np.any(numerator_values <= 0)
        or np.any(denominator_values <= 0)
    ):
        raise SpeedParityError("Paired timings must be positive and finite")
    return float(
        np.exp(np.median(np.log(numerator_values / denominator_values)))
    )


def ratio_gate(
    observed: float,
    threshold: float,
    *,
    calibration_mode: bool,
) -> dict[str, bool]:
    """Evaluate an inclusive ratio target without blocking during calibration."""
    within_limit = bool(observed <= threshold)
    return {
        "within_limit": within_limit,
        "blocking_pass": bool(calibration_mode or within_limit),
    }


def summarize_attempt(
    records: list[dict[str, Any]],
    run_name: str,
    params: dict[str, Any],
    speed_params: dict[str, Any],
    speed_cfg: dict[str, Any],
    attempt: int,
    expected_fingerprints: dict[str, str] | None = None,
) -> dict[str, Any]:
    """Aggregate measured paired trials and evaluate noise and ratio policies."""
    frame = validate_trial_records(
        records,
        run_name,
        params,
        speed_params,
        speed_cfg,
        attempt,
        expected_fingerprints,
    )
    measured = frame.loc[frame["phase"] == "measured"].copy()
    measured["cell"] = np.where(
        measured["engine"] == "r_deseq2",
        "r_deseq2",
        "pydeseq2_" + measured["python_input_mode"].astype(str),
    )
    core = measured.pivot(
        index="repetition",
        columns="cell",
        values="core_seconds",
    ).sort_index()
    required_cells = {"r_deseq2", "pydeseq2_explicit"}
    if not required_cells.issubset(core.columns) or core[list(required_cells)].isna().any().any():
        raise SpeedParityError(f"{run_name} lacks complete primary timing pairs")

    cell_summaries: dict[str, dict[str, dict[str, float]]] = {}
    noisy_cells: list[str] = []
    for cell in sorted(core.columns):
        cell_rows = measured.loc[measured["cell"] == cell]
        phase_summaries = {
            timing_field: robust_timing_summary(
                cell_rows[timing_field].to_numpy()
            )
            for timing_field in (
                "construction_seconds",
                "fit_seconds",
                "results_seconds",
                "core_seconds",
                "process_wall_seconds",
            )
        }
        cell_summaries[cell] = phase_summaries
        if (
            phase_summaries["core_seconds"]["robust_cv"]
            > speed_cfg["noise"]["max_robust_cv"]
        ):
            noisy_cells.append(cell)

    primary_ratio = paired_median_ratio(
        core["pydeseq2_explicit"], core["r_deseq2"]
    )
    provisional = ratio_gate(
        primary_ratio,
        speed_params["provisional_ratio_max"],
        calibration_mode=speed_cfg["calibration_mode"],
    )
    aspirational = ratio_gate(
        primary_ratio,
        speed_cfg["aspirational_ratio_max"],
        calibration_mode=True,
    )
    adata_ratio: float | None = None
    adata_gate: dict[str, bool] | None = None
    if speed_params.get("benchmark_adata", False):
        if "pydeseq2_anndata" not in core:
            raise SpeedParityError(f"{run_name} lacks AnnData timing pairs")
        adata_ratio = paired_median_ratio(
            core["pydeseq2_anndata"], core["pydeseq2_explicit"]
        )
        adata_gate = ratio_gate(
            adata_ratio,
            speed_params["provisional_adata_overhead_max"],
            calibration_mode=speed_cfg["calibration_mode"],
        )

    blocking_pass = provisional["blocking_pass"] and (
        adata_gate is None or adata_gate["blocking_pass"]
    )
    return {
        "run_name": run_name,
        "dataset": params["dataset"],
        "attempt": attempt,
        "measured_repetitions": speed_cfg["measured_repetitions"],
        "pydeseq2_over_r_core_ratio": primary_ratio,
        "provisional_ratio_max": speed_params["provisional_ratio_max"],
        "provisional_ratio_within_limit": provisional["within_limit"],
        "aspirational_ratio_max": speed_cfg["aspirational_ratio_max"],
        "aspirational_ratio_within_limit": aspirational["within_limit"],
        "anndata_over_explicit_core_ratio": adata_ratio,
        "provisional_adata_overhead_max": speed_params.get(
            "provisional_adata_overhead_max"
        ),
        "provisional_adata_within_limit": (
            None if adata_gate is None else adata_gate["within_limit"]
        ),
        "calibration_mode": speed_cfg["calibration_mode"],
        "blocking_ratio_pass": blocking_pass,
        "noisy": bool(noisy_cells),
        "noisy_cells": ",".join(noisy_cells),
        "cell_summaries": cell_summaries,
    }


def should_retry_noisy(
    summary: dict[str, Any],
    attempt: int,
    speed_cfg: dict[str, Any],
) -> bool:
    """Return whether the complete dataset should receive its one noise retry."""
    return bool(
        summary["noisy"]
        and speed_cfg["noise"]["retry_noisy_once"]
        and attempt == 1
    )


def final_run_status(
    summary: dict[str, Any],
    speed_cfg: dict[str, Any],
) -> str:
    """Classify valid evidence without calling infrastructure noise a regression."""
    if summary["noisy"]:
        return "inconclusive_infrastructure"
    if not summary["blocking_ratio_pass"]:
        return "failed"
    return (
        "calibration_observation"
        if speed_cfg["calibration_mode"]
        else "passed"
    )


def validate_parity_precondition(
    parity_cfg: dict[str, Any],
    run_names: list[str],
    pydeseq2_source: dict[str, Any],
) -> dict[str, Any]:
    """Require fresh, passing numerical parity evidence before benchmarking."""
    output_root = run_parity.repo_path(parity_cfg["output_dir"])
    summary_path = output_root / "parity_summary.tsv"
    if not summary_path.is_file():
        raise SpeedParityError("Numerical parity summary is missing")
    summary = pd.read_csv(summary_path, sep="\t", dtype=str, keep_default_na=False)
    if summary["run_name"].duplicated().any():
        raise SpeedParityError("Numerical parity summary contains duplicate runs")
    summary = summary.set_index("run_name")
    config_sha = run_parity.sha256_file(run_parity.CONFIG_PATH)
    current_versions = validate_python_environment(pydeseq2_source, parity_cfg)
    evidence: dict[str, Any] = {}
    for run_name in run_names:
        if run_name not in summary.index:
            raise SpeedParityError(f"Numerical parity result is missing {run_name}")
        row = summary.loc[run_name]
        if row["status"] != "passed" or row["passed"].lower() != "true":
            raise SpeedParityError(f"Numerical parity did not pass for {run_name}")
        run_dir = output_root / run_name
        gates = pd.read_csv(run_dir / "gate_results.tsv", sep="\t")
        if gates.empty:
            raise SpeedParityError(
                f"Numerical parity gate evidence is empty for {run_name}"
            )
        if gates["gate"].duplicated().any():
            raise SpeedParityError(
                f"Numerical parity gate evidence is duplicated for {run_name}"
            )
        try:
            expected_gate_count = int(row["gates"])
        except (KeyError, TypeError, ValueError) as error:
            raise SpeedParityError(
                f"Numerical parity gate count is invalid for {run_name}"
            ) from error
        if len(gates) != expected_gate_count:
            raise SpeedParityError(
                f"Numerical parity gate count changed for {run_name}"
            )
        passed = gates["passed"].astype(str).str.lower().eq("true")
        if not passed.all():
            raise SpeedParityError(f"Numerical parity gates failed for {run_name}")
        provenance = json.loads(
            (run_dir / "provenance.json").read_text(encoding="utf-8")
        )
        if provenance.get("config_sha256") != config_sha:
            raise SpeedParityError(f"Numerical parity config is stale for {run_name}")
        observed_source = provenance.get("pydeseq2", {})
        if (
            observed_source.get("git_commit") != pydeseq2_source["git_commit"]
            or observed_source.get("module_path") != pydeseq2_source["module_path"]
            or observed_source.get("package_version")
            != pydeseq2_source["package_version"]
        ):
            raise SpeedParityError(
                f"Numerical parity source is stale for {run_name}"
            )
        if provenance.get("observed_python_versions") != current_versions:
            raise SpeedParityError(
                f"Numerical parity Python environment is stale for {run_name}"
            )
        expected_inputs = parity_cfg["task_runs"][run_name].get(
            "input_hashes", {}
        )
        if provenance.get("inputs") != expected_inputs:
            raise SpeedParityError(
                f"Numerical parity input hashes changed for {run_name}"
            )
        mode = parity_cfg["task_runs"][run_name]["mode"]
        fingerprints = {
            "r_deseq2": fingerprint_outputs(
                _scientific_output_paths(run_dir, "r", mode)
            ),
            "pydeseq2_explicit": fingerprint_outputs(
                _scientific_output_paths(run_dir, "py", mode)
            ),
        }
        if mode == "tximport":
            fingerprints["pydeseq2_anndata"] = fingerprint_outputs(
                _scientific_output_paths(run_dir, "py_pytximport", mode)
            )
        evidence[run_name] = {
            "provenance": provenance,
            "scientific_fingerprints": fingerprints,
        }
    return evidence


def _host_provenance() -> dict[str, Any]:
    cpu_model = ""
    cpuinfo = Path("/proc/cpuinfo")
    if cpuinfo.is_file():
        for line in cpuinfo.read_text(encoding="utf-8").splitlines():
            if line.startswith("model name"):
                cpu_model = line.split(":", 1)[1].strip()
                break
    affinity = (
        sorted(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else []
    )
    try:
        blas = np.show_config(mode="dicts").get("Build Dependencies", {}).get(
            "blas", {}
        )
    except TypeError:
        blas = {}
    return {
        "platform": platform.platform(),
        "machine": platform.machine(),
        "python": platform.python_version(),
        "cpu_model": cpu_model,
        "logical_cpus": os.cpu_count(),
        "cpu_affinity": affinity,
        "load_average": list(os.getloadavg()) if hasattr(os, "getloadavg") else [],
        "numpy_blas": blas,
        "thread_environment": validate_thread_environment(
            SPEED_CFG["native_math_threads"]
        ),
    }


def _flatten_summary(summary: dict[str, Any]) -> dict[str, Any]:
    row = {key: value for key, value in summary.items() if key != "cell_summaries"}
    for cell, phases in summary["cell_summaries"].items():
        for phase, metrics in phases.items():
            for metric, value in metrics.items():
                row[f"{cell}_{phase}_{metric}"] = value
    return row


def configure_logging(output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=False)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(output_dir / "run_speed_parity.log"),
            logging.StreamHandler(sys.stdout),
        ],
        force=True,
    )


def run_speed_suite(
    speed_cfg: dict[str, Any],
) -> tuple[int, Path]:
    """Run every enabled dataset and write one immutable invocation bundle."""
    validate_speed_config(speed_cfg)
    validate_thread_environment(speed_cfg["native_math_threads"])
    checkout = os.environ.get("PYDESEQ2_REPO")
    if not checkout:
        raise SpeedParityError("PYDESEQ2_REPO must point to the checkout under test")
    pydeseq2_source = run_parity.verify_pydeseq2_checkout(Path(checkout))
    python_versions = validate_python_environment(pydeseq2_source, PARITY_CFG)
    configured = configured_speed_runs({**PARITY_CFG, "speed": speed_cfg})
    enabled = [name for name, item in configured.items() if item["speed"]["run"]]

    invocation_id = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    output_root = run_parity.repo_path(speed_cfg["output_dir"])
    invocation_dir = output_root / invocation_id
    configure_logging(invocation_dir)
    logger = logging.getLogger(__name__)
    logger.info("Benchmarking PyDESeq2 commit %s", pydeseq2_source["git_commit"])

    parity_evidence = validate_parity_precondition(
        PARITY_CFG, enabled, pydeseq2_source
    )
    all_records: list[dict[str, Any]] = []
    attempt_rows: list[dict[str, Any]] = []
    run_results: list[dict[str, Any]] = []
    failures = 0

    def persist_record(record: dict[str, Any]) -> None:
        all_records.append(record)
        pd.DataFrame(all_records, columns=RAW_COLUMNS).to_csv(
            invocation_dir / "speed_trials.tsv", sep="\t", index=False
        )

    for run_name, item in configured.items():
        params = item["analysis"]
        speed_params = item["speed"]
        if not speed_params["run"]:
            logger.info("Skipping %s because run=false", run_name)
            run_results.append(
                {
                    "run_name": run_name,
                    "dataset": params["dataset"],
                    "status": "skipped",
                    "selected_attempt": None,
                }
            )
            continue
        logger.info("Benchmarking %s", run_name)
        try:
            selected: dict[str, Any] | None = None
            max_attempts = (
                2 if speed_cfg["noise"]["retry_noisy_once"] else 1
            )
            for attempt in range(1, max_attempts + 1):
                records = run_attempt(
                    run_name,
                    params,
                    speed_params,
                    speed_cfg,
                    attempt,
                    invocation_dir / run_name,
                    record_callback=persist_record,
                    expected_source=pydeseq2_source,
                )
                summary = summarize_attempt(
                    records,
                    run_name,
                    params,
                    speed_params,
                    speed_cfg,
                    attempt,
                    parity_evidence[run_name]["scientific_fingerprints"],
                )
                attempt_rows.append(summary)
                selected = summary
                if not should_retry_noisy(summary, attempt, speed_cfg):
                    break
                logger.warning(
                    "%s attempt %d is noisy in %s",
                    run_name,
                    attempt,
                    summary["noisy_cells"],
                )
            assert selected is not None
            status = final_run_status(selected, speed_cfg)
            if status in {"inconclusive_infrastructure", "failed"}:
                failures += 1
            run_results.append(
                {
                    "run_name": run_name,
                    "dataset": params["dataset"],
                    "status": status,
                    "selected_attempt": selected["attempt"],
                    **{
                        key: value
                        for key, value in selected.items()
                        if key != "cell_summaries"
                    },
                }
            )
        except Exception as error:
            logger.exception("%s speed benchmark failed", run_name)
            failures += 1
            run_results.append(
                {
                    "run_name": run_name,
                    "dataset": params["dataset"],
                    "status": "error",
                    "selected_attempt": None,
                    "error": str(error),
                }
            )

        pd.DataFrame(all_records, columns=RAW_COLUMNS).to_csv(
            invocation_dir / "speed_trials.tsv", sep="\t", index=False
        )
        pd.DataFrame([_flatten_summary(row) for row in attempt_rows]).to_csv(
            invocation_dir / "speed_attempts.tsv", sep="\t", index=False
        )
        pd.DataFrame(run_results).to_csv(
            invocation_dir / "speed_summary.tsv", sep="\t", index=False
        )

    report = {
        "invocation_id": invocation_id,
        "status": "failed" if failures else "passed",
        "speed_config": speed_cfg,
        "runs": run_results,
        "attempts": attempt_rows,
    }
    (invocation_dir / "speed_summary.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    provenance = {
        "invocation_id": invocation_id,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "command": sys.argv,
        "config_path": str(run_parity.CONFIG_PATH),
        "config_sha256": run_parity.sha256_file(run_parity.CONFIG_PATH),
        "pydeseq2": pydeseq2_source,
        "python_versions": python_versions,
        "host": _host_provenance(),
        "parity_precondition": parity_evidence,
    }
    (invocation_dir / "provenance.json").write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    output_root.mkdir(parents=True, exist_ok=True)
    (output_root / "latest_run.json").write_text(
        json.dumps(
            {
                "invocation_id": invocation_id,
                "path": str(invocation_dir),
                "status": report["status"],
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    return (1 if failures else 0), invocation_dir


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--run-name", help=argparse.SUPPRESS)
    parser.add_argument(
        "--input-mode",
        choices=("explicit", "anndata"),
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--worker-output-dir",
        type=Path,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--warmup-repetitions",
        type=int,
        help="Override the configured warm-up count for this invocation.",
    )
    parser.add_argument(
        "--measured-repetitions",
        type=int,
        help="Override the configured measured count for this invocation.",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if args.worker:
        if not args.run_name or not args.input_mode or args.worker_output_dir is None:
            raise SpeedParityError(
                "Worker mode requires --run-name, --input-mode, and "
                "--worker-output-dir"
            )
        return run_python_worker(
            args.run_name, args.input_mode, args.worker_output_dir
        )
    if args.run_name or args.input_mode or args.worker_output_dir is not None:
        raise SpeedParityError("Worker-only arguments require --worker")
    speed_cfg = deepcopy(SPEED_CFG)
    if args.warmup_repetitions is not None:
        speed_cfg["warmup_repetitions"] = args.warmup_repetitions
    if args.measured_repetitions is not None:
        speed_cfg["measured_repetitions"] = args.measured_repetitions
    exit_code, invocation_dir = run_speed_suite(speed_cfg)
    logging.getLogger(__name__).info("Speed evidence: %s", invocation_dir)
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
