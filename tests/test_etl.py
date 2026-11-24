"""Tests: data-integrity checks, cleaning rules, graph shape, stats, CLI."""
import json
import math

import pytest

from bioassay_etl.cli import main
from bioassay_etl.graph import build_graph
from bioassay_etl.ingest import clean_records, cleaning_accuracy
from bioassay_etl.pipeline import run_pipeline
from bioassay_etl.spark_job import local_activity_counts
from bioassay_etl.stats import describe, eda_report, integrity_report, pearson
from bioassay_etl.synthetic import generate_records


def row(**kw):
    base = {"assay_id": "A1", "compound_id": "C1", "target_id": "T1", "cell_line": "HeLa",
            "assay_type": "IC50", "value": 100, "unit": "nM", "relation": "=", "source": "internal"}
    return {**base, **kw}


# --- cleaning -----------------------------------------------------------
def test_unit_conversion_and_pic50():
    clean, q, _ = clean_records([row(value=1, unit="uM"), row(assay_id="A2", value=1, unit="mM")])
    assert not q
    assert clean[0]["value_nm"] == 1000.0 and clean[0]["pic50"] == 6.0
    assert clean[1]["value_nm"] == 1_000_000.0


def test_each_defect_is_quarantined_with_reason():
    records = [row(), row(assay_id="A2", compound_id=""), row(assay_id="A3", value=-5),
               row(assay_id="A4", unit="furlongs"), row(assay_id="A5", assay_type="XYZ"),
               row(assay_id="A6")]  # duplicate measurement of row()
    clean, q, summary = clean_records(records)
    assert len(clean) == 1 and summary["quarantined"] == 5
    assert {r["quarantine_reason"] for r in q} == {
        "missing_id", "bad_value", "bad_unit", "bad_assay_type", "duplicate"}


def test_out_of_range_quarantined():
    clean, q, _ = clean_records([row(value=1e12, unit="nM")])
    assert not clean and q[0]["quarantine_reason"] == "out_of_range"


def test_cleaning_accuracy_helper():
    assert cleaning_accuracy(1000, 10, 10) == 1.0
    assert cleaning_accuracy(1000, 0, 0) == 1.0
    assert cleaning_accuracy(1000, 9, 10) == pytest.approx(0.9)


# --- graph ---------------------------------------------------------------
def test_graph_nodes_edges_and_activity():
    records, _ = generate_records(500, seed=7)
    clean, _, _ = clean_records(records)
    g = build_graph(clean)
    s = g.stats()
    assert s["nodes_by_type"]["Compound"] > 0 and s["nodes_by_type"]["Target"] == 40
    assert s["edges_by_type"]["MEASURES"] == len(clean)
    for cid, tid in [(e["src"], e["dst"]) for e in g.edges if e["type"] == "ACTIVE_AGAINST"]:
        assert cid.startswith("compound:") and tid.startswith("target:")
        assert cid in g.neighbors(tid)
    assert g.top_hubs(3, "Target")


def test_graph_export(tmp_path):
    clean, _, _ = clean_records([row()])
    g = build_graph(clean)
    paths = g.export(tmp_path)
    data = json.loads(open(paths["graph_json"], encoding="utf-8").read())
    assert data["nodes"] and data["edges"]


# --- stats / integrity ----------------------------------------------------
def test_stats_and_integrity():
    clean, q, summary = clean_records([row(value=v) for v in (10, 100, 1000, 10000)])
    # distinct values -> no duplicates; relabel assay ids not needed (key excludes assay_id)
    rep = integrity_report(clean, len(q), summary["input"])
    assert rep["passed"] and rep["completeness"] == 1.0
    assert describe([1, 2, 3, 4])["median"] == 2.5
    assert pearson([1, 2, 3], [2, 4, 6]) == 1.0
    eda = eda_report(clean)
    assert eda["corr_log10_value_vs_pic50"] == -1.0  # by construction, sanity check
    assert 0.0 <= eda["active_rate"] <= 1.0


def test_integrity_fails_on_row_mismatch():
    clean, _, _ = clean_records([row()])
    assert not integrity_report(clean, 5, 99)["passed"]


def test_local_activity_counts_matches_graph_inputs():
    clean, _, _ = clean_records([row(), row(assay_id="A2", value=100000)])
    counts = local_activity_counts(clean)
    assert counts == [{"compound_id": "C1", "target_id": "T1", "measurements": 2, "actives": 1}]


# --- pipeline / CLI --------------------------------------------------------
def test_pipeline_end_to_end(tmp_path):
    report = run_pipeline(None, tmp_path / "out", records=[row(), row(assay_id="A2", value=-1)])
    assert report["integrity"]["passed"] and report["cleaning"]["clean"] == 1
    assert (tmp_path / "out" / "report.json").exists()
    assert (tmp_path / "out" / "graph" / "nodes.csv").exists()


def test_cli_bench_and_demo(tmp_path, capsys):
    assert main(["bench", "--rows", "2000"]) == 0
    assert "cleaning_accuracy=" in capsys.readouterr().out
    assert main(["demo", "--rows", "500", "--out", str(tmp_path / "demo")]) == 0


def test_generate_records_rejects_bad_args():
    with pytest.raises(ValueError):
        generate_records(0)
    with pytest.raises(ValueError):
        generate_records(10, defect_rate=0.9)
