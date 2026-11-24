"""Deterministic synthetic bioassay data with known ground truth.

The generator first creates *valid* records, then (optionally) injects a
known number of defects (duplicates, bad units, missing ids, non-positive
values). Because the injected defect count is known, the benchmark can
score cleaning accuracy exactly instead of estimating it.
"""
from __future__ import annotations

import random
from typing import Any

ASSAY_TYPES = ["IC50", "EC50", "Ki", "Kd"]
TARGETS = [f"TGT{i:03d}" for i in range(1, 41)]
CELL_LINES = ["HEK293", "HeLa", "A549", "MCF7", "HepG2", "CHO-K1", None]
UNITS = ["nM", "uM", "mM"]  # weighted by position in _pick_unit


def _pick_unit(rng: random.Random) -> str:
    r = rng.random()
    return "nM" if r < 0.7 else ("uM" if r < 0.95 else "mM")


def generate_records(n: int, seed: int = 42, defect_rate: float = 0.0) -> tuple[list[dict[str, Any]], dict[str, int]]:
    """Return (records, ground_truth) where ground_truth counts injected defects."""
    if n <= 0:
        raise ValueError("n must be positive")
    if not 0.0 <= defect_rate < 0.5:
        raise ValueError("defect_rate must be in [0, 0.5)")
    rng = random.Random(seed)
    records: list[dict[str, Any]] = []
    for i in range(n):
        # Log-normal-ish potency in nM, converted to a display unit.
        value_nm = round(10 ** rng.gauss(2.0, 1.1), 4)  # ~100 nM median
        unit = _pick_unit(rng)
        factor = {"nM": 1.0, "uM": 1e3, "mM": 1e6}[unit]
        records.append({
            "assay_id": f"ASY{i:07d}",
            "compound_id": f"CMP{rng.randint(1, max(10, n // 5)):05d}",
            "target_id": rng.choice(TARGETS),
            "cell_line": rng.choice(CELL_LINES),
            "assay_type": rng.choice(ASSAY_TYPES),
            "value": round(value_nm / factor, 10),  # 10 dp: 1 nM = 1e-6 mM must not round to 0
            "unit": unit,
            "relation": rng.choice(["=", "=", "=", "<", ">"]),
            "source": rng.choice(["internal", "public", "partner"]),
        })

    truth = {"injected_defects": 0, "valid_records": n}
    if defect_rate:
        k = int(n * defect_rate)
        kinds = ["duplicate", "bad_unit", "missing_id", "non_positive"]
        planned = [kinds[j % len(kinds)] for j in range(k)]
        mutated: set[int] = set()
        # Phase 1: in-place mutations (each on a distinct row).
        for kind in planned:
            if kind == "duplicate":
                continue
            idx = rng.randrange(len(records))
            while idx in mutated:
                idx = rng.randrange(len(records))
            mutated.add(idx)
            if kind == "bad_unit":
                records[idx] = {**records[idx], "unit": "furlongs"}
            elif kind == "missing_id":
                records[idx] = {**records[idx], "compound_id": ""}
            else:
                records[idx] = {**records[idx], "value": -abs(records[idx]["value"]) or -1.0}
            truth["injected_defects"] += 1
        # Phase 2: duplicates copy only still-valid rows, so a duplicated
        # row can never be "rescued" by later mutating the original.
        for kind in planned:
            if kind != "duplicate":
                continue
            idx = rng.randrange(len(records))
            while idx in mutated:
                idx = rng.randrange(len(records))
            records.append(dict(records[idx]))
            truth["injected_defects"] += 1
    return records, truth
