"""Command line interface.

  python -m bioassay_etl.cli generate --rows 10000 --out data/raw.csv
  python -m bioassay_etl.cli run --input data/raw.csv --out out/
  python -m bioassay_etl.cli bench --rows 200000
  python -m bioassay_etl.cli demo --rows 5000 --out out/
"""
from __future__ import annotations

import argparse
import json
import sys
import time

from .ingest import cleaning_accuracy, clean_records, write_csv
from .pipeline import run_pipeline
from .synthetic import generate_records
from .graph import build_graph


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="bioassay-etl", description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)

    g = sub.add_parser("generate", help="write a synthetic raw CSV")
    g.add_argument("--rows", type=int, default=10_000)
    g.add_argument("--seed", type=int, default=42)
    g.add_argument("--defect-rate", type=float, default=0.002)
    g.add_argument("--out", default="data/raw.csv")

    r = sub.add_parser("run", help="run the ETL pipeline on a CSV/JSON file")
    r.add_argument("--input", required=True)
    r.add_argument("--out", default="out")

    b = sub.add_parser("bench", help="benchmark on synthetic data with known defects")
    b.add_argument("--rows", type=int, default=200_000)
    b.add_argument("--seed", type=int, default=42)
    b.add_argument("--defect-rate", type=float, default=0.002)

    d = sub.add_parser("demo", help="generate synthetic data and run the full pipeline")
    d.add_argument("--rows", type=int, default=5_000)
    d.add_argument("--seed", type=int, default=42)
    d.add_argument("--out", default="out")

    args = ap.parse_args(argv)

    if args.cmd == "generate":
        records, truth = generate_records(args.rows, args.seed, args.defect_rate)
        write_csv(records, args.out, ["assay_id", "compound_id", "target_id", "cell_line",
                                      "assay_type", "value", "unit", "relation", "source"])
        print(f"wrote {len(records)} rows ({truth['injected_defects']} injected defects) -> {args.out}")
        return 0

    if args.cmd == "run":
        report = run_pipeline(args.input, args.out)
        print(json.dumps({k: report[k] for k in ("cleaning", "integrity", "graph", "elapsed_s")}, indent=1))
        return 0 if report["integrity"]["passed"] else 1

    if args.cmd == "demo":
        records, truth = generate_records(args.rows, args.seed, defect_rate=0.002)
        report = run_pipeline(None, args.out, records=records)
        print(json.dumps({k: report[k] for k in ("cleaning", "integrity", "graph", "elapsed_s")}, indent=1))
        return 0 if report["integrity"]["passed"] else 1

    if args.cmd == "bench":
        t0 = time.perf_counter()
        records, truth = generate_records(args.rows, args.seed, args.defect_rate)
        t_gen = time.perf_counter()
        clean, quarantine, summary = clean_records(records)
        t_clean = time.perf_counter()
        graph = build_graph(clean)
        t_graph = time.perf_counter()
        acc = cleaning_accuracy(summary["input"], len(quarantine), truth["injected_defects"])
        total = t_graph - t0
        print(f"rows_input={summary['input']} clean={len(clean)} quarantined={len(quarantine)} "
              f"injected_defects={truth['injected_defects']}")
        print(f"cleaning_accuracy={acc:.4%} nodes={len(graph.nodes)} edges={len(graph.edges)}")
        print(f"generate={t_gen - t0:.2f}s clean={t_clean - t_gen:.2f}s graph={t_graph - t_clean:.2f}s total={total:.2f}s")
        print(f"throughput={summary['input'] / total:,.0f} rows/s")
        return 0 if acc >= 0.999 else 1
    return 2


if __name__ == "__main__":
    sys.exit(main())
