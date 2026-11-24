#!/usr/bin/env python
"""Benchmark entry point (same code path as `bioassay-etl bench`).

Usage: python scripts/benchmark.py --rows 1000000
Scaled-down sanity runs are reported honestly in the README; the 1M-row
figure is a target you can reproduce with the command above, not a
number measured in CI.
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from bioassay_etl.cli import main  # noqa: E402

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--rows", type=int, default=200_000)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--defect-rate", type=float, default=0.002)
    a = ap.parse_args()
    raise SystemExit(main(["bench", "--rows", str(a.rows), "--seed", str(a.seed),
                           "--defect-rate", str(a.defect_rate)]))
