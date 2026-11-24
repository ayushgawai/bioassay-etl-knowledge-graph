"""Pipeline orchestration: ingest -> clean -> graph -> EDA -> report."""
from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any

from .graph import build_graph
from .ingest import clean_records, load_records, write_csv
from .stats import eda_report, integrity_report

CLEAN_FIELDS = ["assay_id", "compound_id", "target_id", "cell_line", "assay_type",
                "value_nm", "pic50", "relation", "source", "is_active"]


def run_pipeline(input_path: str | Path | None, out_dir: str | Path,
                 records: list[dict[str, Any]] | None = None) -> dict[str, Any]:
    t0 = time.perf_counter()
    raw = records if records is not None else load_records(input_path)  # type: ignore[arg-type]
    clean, quarantine, summary = clean_records(raw)
    graph = build_graph(clean)
    out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)
    write_csv(clean, out / "clean.csv", CLEAN_FIELDS)
    if quarantine:
        write_csv(quarantine, out / "quarantine.csv")
    graph_paths = graph.export(out / "graph")
    report = {
        "cleaning": summary,
        "integrity": integrity_report(clean, len(quarantine), len(raw)),
        "eda": eda_report(clean),
        "graph": graph.stats(),
        "outputs": {"clean_csv": str(out / "clean.csv"), "quarantine_csv": str(out / "quarantine.csv"),
                    **graph_paths},
        "elapsed_s": round(time.perf_counter() - t0, 3),
    }
    (out / "report.json").write_text(json.dumps(report, indent=1), encoding="utf-8")
    return report
