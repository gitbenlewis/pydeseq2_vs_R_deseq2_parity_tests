from __future__ import annotations

import importlib.util
import json
from copy import deepcopy
from pathlib import Path

import numpy as np
import pandas as pd
import pytest


SCRIPT_PATH = (
    Path(__file__).resolve().parents[1] / "scripts" / "run_speed_parity.py"
)
SPEC = importlib.util.spec_from_file_location("run_speed_parity", SCRIPT_PATH)
assert SPEC is not None and SPEC.loader is not None
run_speed = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(run_speed)


def _speed_config(
    *,
    run_name: str = "pasilla",
    warmups: int = 1,
    measured: int = 3,
    benchmark_adata: bool = False,
) -> tuple[dict, dict, dict]:
    config = deepcopy(run_speed.SPEED_CFG)
    config["warmup_repetitions"] = warmups
    config["measured_repetitions"] = measured
    config["task_runs"] = {
        run_name: deepcopy(config["task_runs"][run_name])
    }
    config["task_runs"][run_name]["benchmark_adata"] = benchmark_adata
    if benchmark_adata:
        config["task_runs"][run_name][
            "provisional_adata_overhead_max"
        ] = 1.15
    analysis = run_speed.run_parity.configured_runs(
        run_speed.PARITY_CFG
    )[run_name]
    speed_params = run_speed.configured_speed_runs(
        {**run_speed.PARITY_CFG, "speed": config}
    )[run_name]["speed"]
    return config, analysis, speed_params


def _complete_records(
    run_name: str,
    analysis: dict,
    speed_params: dict,
    speed_cfg: dict,
    *,
    primary_ratio: float = 2.0,
    adata_ratio: float = 1.05,
    attempt: int = 1,
) -> list[dict]:
    records = []
    schedule = run_speed.build_schedule(
        run_name,
        speed_cfg["warmup_repetitions"],
        speed_cfg["measured_repetitions"],
        attempt,
        speed_params.get("benchmark_adata", False),
        speed_cfg["random_seed"],
    )
    for item in schedule:
        r_core = 1.0 + 0.1 * item["repetition"]
        if item["engine"] == "r_deseq2":
            core = r_core
            fingerprint = "r-fingerprint"
        elif item["python_input_mode"] == "explicit":
            core = r_core * primary_ratio
            fingerprint = "py-explicit-fingerprint"
        else:
            core = r_core * primary_ratio * adata_ratio
            fingerprint = "py-anndata-fingerprint"
        record = {
            "run_name": run_name,
            "dataset": analysis["dataset"],
            "engine": item["engine"],
            "python_input_mode": item["python_input_mode"],
            "construction_seconds": 0.1,
            "fit_seconds": core - 0.2,
            "results_seconds": 0.2,
            "core_seconds": core,
            "process_wall_seconds": core + 1.0,
            "genes": analysis["expected_genes"],
            "samples": analysis["expected_samples"],
            "dispersion_fit_type": analysis["fit_type"],
            "fingerprint": fingerprint,
            "worker_output_dir": "/tmp/worker",
            **item,
        }
        records.append(record)
    return records


def _write_r_metadata(path: Path, params: dict, **overrides) -> None:
    values = {
        "dataset": params["dataset"],
        "mode": params["mode"],
        "genes": params["expected_genes"],
        "samples": params["expected_samples"],
        "design": params["design"],
        "contrast": (
            f"{params['contrast_factor']}_{params['contrast_numerator']}"
            f"_vs_{params['contrast_denominator']}"
        ),
        "dispersion_fit_type": params["fit_type"],
        "size_factor_type": params["size_factor_fit_type"],
        "n_cpus": 1,
        "parallel": "FALSE",
        "construction_seconds": 0.2,
        "fit_seconds": 1.0,
        "results_seconds": 0.3,
    }
    values.update(overrides)
    pd.DataFrame(values.items(), columns=["key", "value"]).to_csv(
        path, sep="\t", index=False
    )


def test_speed_config_references_known_parity_runs_and_honors_run_false():
    parity = deepcopy(run_speed.PARITY_CFG)
    configured = run_speed.configured_speed_runs(parity)

    assert set(configured) == set(parity["speed"]["task_runs"])
    assert set(configured).issubset(parity["task_runs"])

    parity["task_runs"]["pasilla"]["run"] = False
    assert (
        run_speed.configured_speed_runs(parity)["pasilla"]["speed"]["run"]
        is False
    )

    parity = deepcopy(run_speed.PARITY_CFG)
    parity["speed"]["task_runs"]["pasilla"]["run"] = False
    assert (
        run_speed.configured_speed_runs(parity)["pasilla"]["speed"]["run"]
        is False
    )


@pytest.mark.parametrize(
    ("path", "value", "message"),
    [
        (("measured_repetitions",), 0, "measured_repetitions"),
        (("warmup_repetitions",), -1, "warmup_repetitions"),
        (("execution_order",), "fixed", "counterbalanced"),
        (("n_cpus",), 2, "n_cpus"),
        (("native_math_threads",), 2, "native_math_threads"),
        (("aspirational_ratio_max",), np.nan, "aspirational_ratio_max"),
        (("noise", "max_robust_cv"), 1.0, "less than one"),
        (
            ("task_runs", "pasilla", "provisional_ratio_max"),
            0,
            "provisional_ratio_max",
        ),
    ],
)
def test_validate_speed_config_rejects_invalid_execution_settings(
    path, value, message
):
    config = deepcopy(run_speed.SPEED_CFG)
    target = config
    for key in path[:-1]:
        target = target[key]
    target[path[-1]] = value

    with pytest.raises(run_speed.SpeedParityError, match=message):
        run_speed.validate_speed_config(config)


def test_counterbalanced_schedule_is_deterministic_and_excludes_warmups():
    first = run_speed.build_schedule("pasilla", 1, 5, 1, False, 42)
    second = run_speed.build_schedule("pasilla", 1, 5, 1, False, 42)

    assert first == second
    measured = [row for row in first if row["phase"] == "measured"]
    assert len(measured) == 10
    for repetition in range(1, 6):
        pair = [row for row in measured if row["repetition"] == repetition]
        assert {row["engine"] for row in pair} == {
            "r_deseq2",
            "pydeseq2",
        }
    first_engines = [
        next(
            row["engine"]
            for row in measured
            if row["repetition"] == repetition
            and row["order_position"] == 1
        )
        for repetition in range(1, 6)
    ]
    assert all(
        left != right for left, right in zip(first_engines, first_engines[1:])
    )
    assert len([row for row in first if row["phase"] == "warmup"]) == 2


def test_core_seconds_excludes_construction_and_process_wall():
    record = run_speed.timing_record(
        {
            "construction_seconds": 10,
            "fit_seconds": 2,
            "results_seconds": 3,
            "genes": 5,
            "samples": 2,
            "dispersion_fit_type": "parametric",
        },
        run_name="x",
        dataset="x",
        engine="pydeseq2",
        input_mode="explicit",
        process_wall_seconds=30,
        fingerprint="abc",
        worker_output_dir=Path("/tmp/x"),
    )

    assert record["core_seconds"] == 5
    assert record["construction_seconds"] == 10
    assert record["process_wall_seconds"] == 30


def test_scientific_fingerprint_is_prefix_independent(tmp_path):
    left = []
    right = []
    for position in range(4):
        left_path = tmp_path / f"py_{position}.tsv"
        right_path = tmp_path / f"py_pytximport_{position}.tsv"
        left_path.write_text(f"value-{position}\n", encoding="utf-8")
        right_path.write_text(f"value-{position}\n", encoding="utf-8")
        left.append(left_path)
        right.append(right_path)

    assert run_speed.fingerprint_outputs(left) == run_speed.fingerprint_outputs(
        right
    )


def test_r_metadata_parser_validates_identity_and_timings(tmp_path):
    _, params, _ = _speed_config()
    path = tmp_path / "r_metadata.tsv"
    _write_r_metadata(path, params)

    parsed = run_speed.parse_r_metadata(path, params)

    assert parsed["fit_seconds"] == 1.0
    assert parsed["genes"] == params["expected_genes"]

    _write_r_metadata(path, params, fit_seconds="NaN")
    with pytest.raises(run_speed.SpeedParityError, match="finite"):
        run_speed.parse_r_metadata(path, params)

    _write_r_metadata(path, params, dataset="wrong")
    with pytest.raises(run_speed.SpeedParityError, match="dataset"):
        run_speed.parse_r_metadata(path, params)


def test_r_metadata_parser_rejects_duplicate_keys(tmp_path):
    _, params, _ = _speed_config()
    path = tmp_path / "r_metadata.tsv"
    _write_r_metadata(path, params)
    frame = pd.read_csv(path, sep="\t")
    frame = pd.concat([frame, frame.loc[frame["key"] == "fit_seconds"]])
    frame.to_csv(path, sep="\t", index=False)

    with pytest.raises(run_speed.SpeedParityError, match="Duplicate"):
        run_speed.parse_r_metadata(path, params)


def test_trial_validation_rejects_missing_duplicate_and_changed_fingerprint():
    config, params, speed_params = _speed_config()
    records = _complete_records("pasilla", params, speed_params, config)

    run_speed.validate_trial_records(
        records, "pasilla", params, speed_params, config, 1
    )

    with pytest.raises(run_speed.SpeedParityError, match="incomplete"):
        run_speed.validate_trial_records(
            records[:-1], "pasilla", params, speed_params, config, 1
        )

    with pytest.raises(run_speed.SpeedParityError, match="duplicate"):
        run_speed.validate_trial_records(
            records + [records[-1]],
            "pasilla",
            params,
            speed_params,
            config,
            1,
        )

    changed = deepcopy(records)
    changed[-1]["fingerprint"] = "different"
    with pytest.raises(run_speed.SpeedParityError, match="nondeterministic"):
        run_speed.validate_trial_records(
            changed, "pasilla", params, speed_params, config, 1
        )

    expected = {
        "r_deseq2": "r-fingerprint",
        "pydeseq2_explicit": "py-explicit-fingerprint",
    }
    drifted = deepcopy(records)
    for row in drifted:
        if row["engine"] == "pydeseq2":
            row["fingerprint"] = "stable-but-wrong"
    with pytest.raises(run_speed.SpeedParityError, match="numerical parity"):
        run_speed.validate_trial_records(
            drifted,
            "pasilla",
            params,
            speed_params,
            config,
            1,
            expected,
        )


def test_summary_uses_median_log_ratio_and_robust_statistics():
    numerator = pd.Series([1.0, 9.0], index=[1, 2])
    denominator = pd.Series([1.0, 1.0], index=[1, 2])

    assert run_speed.paired_median_ratio(numerator, denominator) == pytest.approx(
        3.0
    )
    summary = run_speed.robust_timing_summary([1.0, 2.0, 3.0])
    assert summary["median"] == 2.0
    assert summary["mad"] == 1.0
    assert summary["robust_cv"] == pytest.approx(1.4826 / 2.0)
    assert summary["iqr"] == 1.0


def test_ratio_gate_is_inclusive_and_calibration_is_nonblocking():
    boundary = run_speed.ratio_gate(
        2.0, 2.0, calibration_mode=False
    )
    above = run_speed.ratio_gate(
        np.nextafter(2.0, np.inf), 2.0, calibration_mode=False
    )
    calibration = run_speed.ratio_gate(
        3.0, 2.0, calibration_mode=True
    )

    assert boundary == {"within_limit": True, "blocking_pass": True}
    assert above == {"within_limit": False, "blocking_pass": False}
    assert calibration == {"within_limit": False, "blocking_pass": True}


def test_noise_retry_and_inconclusive_policy():
    config, _, _ = _speed_config()
    noisy = {"noisy": True, "blocking_ratio_pass": True}
    clean = {"noisy": False, "blocking_ratio_pass": True}

    assert run_speed.should_retry_noisy(noisy, 1, config) is True
    assert run_speed.should_retry_noisy(noisy, 2, config) is False
    assert (
        run_speed.final_run_status(noisy, config)
        == "inconclusive_infrastructure"
    )
    assert run_speed.final_run_status(clean, config) == "calibration_observation"

    config["calibration_mode"] = False
    failed = {"noisy": False, "blocking_ratio_pass": False}
    assert run_speed.final_run_status(failed, config) == "failed"


def test_srp_anndata_is_secondary_to_primary_r_ratio():
    config, params, speed_params = _speed_config(
        run_name="srp254919_tximport",
        benchmark_adata=True,
    )
    records = _complete_records(
        "srp254919_tximport",
        params,
        speed_params,
        config,
        primary_ratio=1.2,
        adata_ratio=1.1,
    )

    summary = run_speed.summarize_attempt(
        records,
        "srp254919_tximport",
        params,
        speed_params,
        config,
        1,
    )

    assert summary["pydeseq2_over_r_core_ratio"] == pytest.approx(1.2)
    assert summary["anndata_over_explicit_core_ratio"] == pytest.approx(1.1)


def test_mocked_attempt_uses_fresh_paths_and_summarizes_without_sleeping():
    config, params, speed_params = _speed_config(measured=2)
    seen_paths = []
    persisted = []

    def trial_runner(run_name, analysis, engine, input_mode, output_dir):
        seen_paths.append(output_dir)
        core = 1.0 if engine == "r_deseq2" else 2.0
        return {
            "run_name": run_name,
            "dataset": analysis["dataset"],
            "engine": engine,
            "python_input_mode": input_mode,
            "construction_seconds": 0.1,
            "fit_seconds": core - 0.2,
            "results_seconds": 0.2,
            "core_seconds": core,
            "process_wall_seconds": core + 1,
            "genes": analysis["expected_genes"],
            "samples": analysis["expected_samples"],
            "dispersion_fit_type": analysis["fit_type"],
            "fingerprint": f"{engine}-{input_mode}",
            "worker_output_dir": str(output_dir),
        }

    records = run_speed.run_attempt(
        "pasilla",
        params,
        speed_params,
        config,
        1,
        Path("/tmp/speed-test"),
        trial_runner=trial_runner,
        record_callback=persisted.append,
    )
    summary = run_speed.summarize_attempt(
        records, "pasilla", params, speed_params, config, 1
    )

    assert len(seen_paths) == len(set(seen_paths)) == 6
    assert persisted == records
    assert sum(row["phase"] == "warmup" for row in records) == 2
    assert summary["measured_repetitions"] == 2
    assert summary["pydeseq2_over_r_core_ratio"] == pytest.approx(2.0)


def test_parity_precondition_rejects_stale_config(tmp_path, monkeypatch):
    config_path = tmp_path / "config.yaml"
    config_path.write_text("current: true\n", encoding="utf-8")
    output_root = tmp_path / "parity"
    run_dir = output_root / "pasilla"
    run_dir.mkdir(parents=True)
    pd.DataFrame(
        [
            {
                "run_name": "pasilla",
                "passed": True,
                "status": "passed",
                "gates": 1,
            }
        ]
    ).to_csv(output_root / "parity_summary.tsv", sep="\t", index=False)
    pd.DataFrame([{"gate": "all", "passed": True}]).to_csv(
        run_dir / "gate_results.tsv", sep="\t", index=False
    )
    for prefix in ("r", "py"):
        for path in run_speed._scientific_output_paths(run_dir, prefix, "matrix"):
            path.write_text("deterministic\n", encoding="utf-8")
    source = {
        "git_commit": "abc123",
        "module_path": "/source/pydeseq2/__init__.py",
        "package_version": run_speed.PARITY_CFG["expected_versions"][
            "pydeseq2"
        ],
    }
    provenance = {
        "config_sha256": run_speed.run_parity.sha256_file(config_path),
        "pydeseq2": source,
        "observed_python_versions": run_speed.run_parity.python_runtime_versions(),
        "inputs": run_speed.PARITY_CFG["task_runs"]["pasilla"][
            "input_hashes"
        ],
    }
    (run_dir / "provenance.json").write_text(
        json.dumps(provenance), encoding="utf-8"
    )
    parity = deepcopy(run_speed.PARITY_CFG)
    parity["output_dir"] = str(output_root)
    monkeypatch.setattr(run_speed.run_parity, "CONFIG_PATH", config_path)

    evidence = run_speed.validate_parity_precondition(
        parity, ["pasilla"], source
    )
    assert (
        evidence["pasilla"]["provenance"]["config_sha256"]
        == provenance["config_sha256"]
    )

    changed_versions = deepcopy(provenance["observed_python_versions"])
    changed_versions["numpy"] = "0.0"
    with monkeypatch.context() as scoped:
        scoped.setattr(
            run_speed.run_parity,
            "python_runtime_versions",
            lambda: changed_versions,
        )
        with pytest.raises(run_speed.SpeedParityError, match="environment is stale"):
            run_speed.validate_parity_precondition(parity, ["pasilla"], source)

    pd.DataFrame(columns=["gate", "passed"]).to_csv(
        run_dir / "gate_results.tsv", sep="\t", index=False
    )
    with pytest.raises(run_speed.SpeedParityError, match="evidence is empty"):
        run_speed.validate_parity_precondition(parity, ["pasilla"], source)
    pd.DataFrame([{"gate": "all", "passed": True}]).to_csv(
        run_dir / "gate_results.tsv", sep="\t", index=False
    )

    provenance["config_sha256"] = "stale"
    (run_dir / "provenance.json").write_text(
        json.dumps(provenance), encoding="utf-8"
    )
    with pytest.raises(run_speed.SpeedParityError, match="config is stale"):
        run_speed.validate_parity_precondition(parity, ["pasilla"], source)


def test_committed_speed_config_is_calibration_only_and_complete():
    config = run_speed.SPEED_CFG

    assert config["calibration_mode"] is True
    assert config["warmup_repetitions"] == 1
    assert config["measured_repetitions"] == 7
    assert set(config["task_runs"]) == {
        "srp254919_tximport",
        "pasilla",
        "pickrell",
    }
    assert all(run["run"] for run in config["task_runs"].values())
    assert config["task_runs"]["srp254919_tximport"][
        "provisional_adata_overhead_max"
    ] == 1.15
