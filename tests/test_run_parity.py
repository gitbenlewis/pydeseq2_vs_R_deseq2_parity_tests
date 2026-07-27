from __future__ import annotations

import hashlib
import importlib.util
from copy import deepcopy
from pathlib import Path

import numpy as np
import pandas as pd
import pytest


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "scripts" / "run_parity.py"
SPEC = importlib.util.spec_from_file_location("run_parity", SCRIPT_PATH)
assert SPEC is not None and SPEC.loader is not None
run_parity = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(run_parity)


def test_merge_params_recursively_applies_run_overrides():
    defaults = {
        "run": True,
        "fit_type": "parametric",
        "gate_overrides": {"result_columns": {"padj": {"max": 1}}},
    }
    overrides = {
        "run": False,
        "gate_overrides": {
            "result_columns": {"padj": {"na_mask_disagreements_max": 1250}}
        },
    }

    merged = run_parity.merge_params(defaults, overrides)

    assert merged["run"] is False
    assert merged["fit_type"] == "parametric"
    assert merged["gate_overrides"]["result_columns"]["padj"] == {
        "max": 1,
        "na_mask_disagreements_max": 1250,
    }


def test_configured_runs_preserves_disabled_run():
    config = {
        "default_params": {"run": True, "n_cpus": 1},
        "task_runs": {
            "enabled": {},
            "disabled": {"run": False},
        },
    }

    runs = run_parity.configured_runs(config)

    assert runs["enabled"] == {"run": True, "n_cpus": 1}
    assert runs["disabled"] == {"run": False, "n_cpus": 1}


def test_validate_run_params_rejects_reversed_contrast():
    params = run_parity.configured_runs(run_parity.PARITY_CFG)["pasilla"]
    params["contrast_numerator"] = "untreated"
    params["contrast_denominator"] = "treated"

    with pytest.raises(run_parity.ParityError, match="test level versus reference"):
        run_parity.validate_run_params(params)


def test_validate_dispersion_fit_type_rejects_fallback():
    with pytest.raises(run_parity.ParityError, match="used dispersion fit"):
        run_parity.validate_dispersion_fit_type(
            observed="mean",
            requested="parametric",
            engine="PyDESeq2",
        )

    run_parity.validate_dispersion_fit_type(
        observed="parametric",
        requested="parametric",
        engine="R DESeq2",
    )


def test_main_does_not_execute_disabled_runs(tmp_path, monkeypatch):
    config = deepcopy(run_parity.PARITY_CFG)
    config["output_dir"] = str(tmp_path / "results")
    for params in config["task_runs"].values():
        params["run"] = False

    monkeypatch.setattr(run_parity, "PARITY_CFG", config)
    monkeypatch.setenv("PYDESEQ2_REPO", str(tmp_path / "source"))
    monkeypatch.setattr(
        run_parity,
        "verify_pydeseq2_checkout",
        lambda _: {
            "git_commit": "abc123",
            "module_path": "/requested/pydeseq2/__init__.py",
            "package_version": config["expected_versions"]["pydeseq2"],
        },
    )

    def unexpected_run(*_args, **_kwargs):
        raise AssertionError("disabled run was executed")

    monkeypatch.setattr(run_parity, "run_named_task", unexpected_run)

    assert run_parity.main() == 0
    summary = pd.read_csv(tmp_path / "results" / "parity_summary.tsv", sep="\t")
    assert summary["status"].tolist() == ["skipped", "skipped", "skipped"]


def test_verify_checksum_rejects_mismatch(tmp_path):
    path = tmp_path / "fixture.tsv"
    path.write_text("gene_id\ts1\ng1\t1\n", encoding="utf-8")

    with pytest.raises(run_parity.ParityError, match="SHA-256 mismatch"):
        run_parity.verify_checksum(path, "0" * 64)

    expected = hashlib.sha256(path.read_bytes()).hexdigest()
    assert run_parity.verify_checksum(path, expected) == expected


def test_load_analysis_inputs_orients_samples_by_genes(tmp_path):
    counts_path = tmp_path / "counts.tsv"
    counts_path.write_text(
        "gene_id\tgene_name\ts2\ts1\n"
        "g1\tone\t2\t1\n"
        "g2\ttwo\t4\t3\n",
        encoding="utf-8",
    )
    lengths_path = tmp_path / "lengths.tsv"
    lengths_path.write_text(
        "gene_id\tgene_name\ts2\ts1\n"
        "g1\tone\t20\t10\n"
        "g2\ttwo\t40\t30\n",
        encoding="utf-8",
    )
    samples_path = tmp_path / "samples.csv"
    samples_path.write_text(
        "sample,treatment\ns1,control\ns2,treated\n",
        encoding="utf-8",
    )

    counts, metadata, lengths = run_parity.load_analysis_inputs(
        counts_path,
        samples_path,
        "sample",
        "treatment",
        ["control", "treated"],
        lengths_path,
    )

    assert counts.index.tolist() == ["s1", "s2"]
    assert counts.columns.tolist() == ["g1", "g2"]
    np.testing.assert_array_equal(counts, [[1, 3], [2, 4]])
    assert metadata["treatment"].cat.categories.tolist() == ["control", "treated"]
    assert lengths is not None
    np.testing.assert_array_equal(lengths, [[10, 30], [20, 40]])


def test_validate_frame_alignment_rejects_reordering():
    left = pd.DataFrame({"s1": [1, 2]}, index=["g1", "g2"])
    right = left.iloc[::-1]

    with pytest.raises(run_parity.ParityError, match="row labels"):
        run_parity.validate_frame_alignment(left, right, "counts")


def test_comparison_metrics_handles_na_and_zero_vs_tiny_pvalue():
    left = pd.DataFrame({"pvalue": [0.0, 0.5, np.nan]}, index=["a", "b", "c"])
    right = pd.DataFrame(
        {"pvalue": [1e-300, 0.5001, np.nan]},
        index=["a", "b", "c"],
    )

    metrics = run_parity.comparison_metrics(left, right)

    assert metrics["finite_pairs"] == 2
    assert metrics["na_mask_disagreements"] == 0
    assert metrics["abs_diff_max"] == pytest.approx(0.0001)
    assert np.isfinite(metrics["spearman"])


def test_comparison_rejects_infinity_instead_of_dropping_it():
    left = pd.DataFrame({"log2FoldChange": [np.inf, 1.0]}, index=["a", "b"])
    right = pd.DataFrame({"log2FoldChange": [-np.inf, 1.0]}, index=["a", "b"])

    with pytest.raises(run_parity.ParityError, match="infinity"):
        run_parity.comparison_metrics(left, right)


def test_na_mask_disagreement_gate_is_inclusive():
    r_outputs = _identical_engine_outputs()
    py_outputs = {
        key: value.copy()
        for key, value in r_outputs.items()
    }
    py_outputs["results"].loc["g2", "padj"] = np.nan
    summary = run_parity.summarize_comparison(r_outputs, py_outputs)
    rules = _boundary_rules()
    rules["result_columns"]["padj"] = {"na_mask_disagreements_max": 1}

    gates = run_parity.evaluate_gates(
        summary,
        r_outputs,
        py_outputs,
        rules,
    )

    padj_gate = next(
        gate
        for gate in gates
        if gate["gate"] == "padj_na_mask_disagreements_max"
    )
    assert padj_gate["observed"] == 1
    assert padj_gate["passed"] is True

    rules["result_columns"]["padj"]["na_mask_disagreements_max"] = 0
    failed = run_parity.evaluate_gates(
        summary,
        r_outputs,
        py_outputs,
        rules,
    )
    assert next(
        gate
        for gate in failed
        if gate["gate"] == "padj_na_mask_disagreements_max"
    )["passed"] is False


def _identical_engine_outputs():
    genes = pd.Index(["g1", "g2", "g3"], name="gene_id")
    samples = pd.Index(["s1", "s2"], name="sample")
    rounded = pd.DataFrame(
        [[1, 2], [3, 4], [5, 6]],
        index=genes,
        columns=samples,
    )
    factors = pd.DataFrame(
        [[0.9, 1.1], [1.0, 1.0], [1.1, 0.9]],
        index=genes,
        columns=samples,
    )
    normalized = rounded / factors
    results = pd.DataFrame(
        {
            "baseMean": normalized.mean(axis=1),
            "log2FoldChange": [-1.0, 0.2, 1.0],
            "lfcSE": [0.2, 0.3, 0.4],
            "stat": [-5.0, 0.67, 2.5],
            "pvalue": [0.001, 0.5, 0.02],
            "padj": [0.003, 0.5, 0.03],
        },
        index=genes,
    )
    return {
        "rounded_counts": rounded,
        "factors": factors,
        "normalized_counts": normalized,
        "results": results,
    }


def _boundary_rules():
    return {
        "factors": {"rtol": 0.0, "atol": 0.0},
        "normalized_counts": {"rtol": 0.0, "atol": 0.0},
        "base_mean": {"rtol": 0.0, "atol": 0.0},
        "result_columns": {
            "log2FoldChange": {
                "pearson_min": 1.0,
                "spearman_min": 1.0,
                "abs_diff_p95_max": 0.0,
                "abs_diff_max": 0.0,
                "na_mask_disagreements_max": 0,
            },
            "pvalue": {
                "pearson_min": 1.0,
                "spearman_min": 1.0,
                "abs_diff_p95_max": 0.0,
                "abs_diff_max": 0.0,
                "na_mask_disagreements_max": 0,
            },
            "padj": {
                "pearson_min": 1.0,
                "spearman_min": 1.0,
                "abs_diff_p95_max": 0.0,
                "abs_diff_max": 0.0,
                "na_mask_disagreements_max": 0,
            },
        },
        "sign_concordance_min": 1.0,
        "significant_jaccard_min": {"0.1": 1.0},
        "significant_sets_exact": [0.05, 0.1],
    }


def test_every_configured_gate_passes_at_its_boundary():
    outputs = _identical_engine_outputs()
    summary = run_parity.summarize_comparison(outputs, outputs)

    gates = run_parity.evaluate_gates(
        summary,
        outputs,
        outputs,
        _boundary_rules(),
    )

    assert gates
    assert all(gate["passed"] for gate in gates)


def test_gate_failure_is_reported_without_hiding_other_gates():
    outputs = _identical_engine_outputs()
    summary = run_parity.summarize_comparison(outputs, outputs)
    rules = _boundary_rules()
    rules["result_columns"]["log2FoldChange"]["pearson_min"] = 1.0001

    gates = run_parity.evaluate_gates(summary, outputs, outputs, rules)

    failed = [gate["gate"] for gate in gates if not gate["passed"]]
    assert failed == ["log2FoldChange_pearson_min"]
    assert len(gates) > 1


def test_allclose_gate_reports_empty_finite_comparison_as_failure():
    frame = pd.DataFrame({"value": [np.nan]}, index=["g1"])
    gates = []

    run_parity._allclose_gate(
        gates,
        "all_na",
        frame,
        frame,
        {"rtol": 0.0, "atol": 0.0},
    )

    assert len(gates) == 1
    assert gates[0]["passed"] is False
    assert np.isnan(gates[0]["observed"])


def test_committed_config_contains_every_approved_gate():
    config = run_parity.PARITY_CFG
    assert list(config["task_runs"]) == [
        "srp254919_tximport",
        "pasilla",
        "pickrell",
    ]
    assert all(params["run"] is True for params in config["task_runs"].values())

    strict = config["gate_profiles"]["tximport_strict"]
    assert strict["result_columns"] == {
        "log2FoldChange": {
            "pearson_min": 0.999999,
            "abs_diff_max": 0.001,
            "na_mask_disagreements_max": 0,
        },
        "lfcSE": {
            "pearson_min": 0.9999,
            "abs_diff_max": 0.02,
        },
        "stat": {
            "pearson_min": 0.9999,
            "abs_diff_max": 0.05,
        },
        "pvalue": {
            "spearman_min": 0.9999,
            "abs_diff_max": 0.005,
            "na_mask_disagreements_max": 0,
        },
        "padj": {
            "spearman_min": 0.9999,
            "abs_diff_max": 0.01,
            "na_mask_disagreements_max": 0,
        },
    }
    assert strict["sign_concordance_min"] == 1.0
    assert strict["significant_sets_exact"] == [0.05, 0.1]

    raw = config["gate_profiles"]["raw_regression"]
    assert raw["result_columns"]["log2FoldChange"] == {
        "pearson_min": 0.9999,
        "spearman_min": 0.9999,
        "abs_diff_p95_max": 5e-4,
        "abs_diff_max": 0.025,
        "na_mask_disagreements_max": 0,
    }
    assert raw["result_columns"]["pvalue"] == {
        "pearson_min": 0.9997,
        "spearman_min": 0.9997,
        "abs_diff_p95_max": 0.02,
        "abs_diff_max": 0.20,
        "na_mask_disagreements_max": 1,
    }
    assert raw["result_columns"]["padj"] == {
        "pearson_min": 0.98,
        "spearman_min": 0.998,
        "abs_diff_p95_max": 0.10,
        "abs_diff_max": 0.25,
        "na_mask_disagreements_max": 0,
    }
    assert raw["sign_concordance_min"] == 0.999
    assert raw["significant_jaccard_min"] == {"0.1": 0.97}
    assert (
        config["task_runs"]["pickrell"]["gate_overrides"]["result_columns"][
            "padj"
        ]["na_mask_disagreements_max"]
        == 1250
    )


def test_pytximport_input_modes_require_numerical_identity():
    explicit = _identical_engine_outputs()
    compatible = {
        key: value.copy()
        for key, value in explicit.items()
    }
    run_parity._assert_python_input_modes_equal(
        explicit,
        compatible,
        rtol=1e-12,
        atol=1e-12,
    )

    compatible["normalized_counts"].iloc[0, 0] += 1e-3
    with pytest.raises(run_parity.ParityError, match="pytximport AnnData"):
        run_parity._assert_python_input_modes_equal(
            explicit,
            compatible,
            rtol=1e-12,
            atol=1e-12,
        )
