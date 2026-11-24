"""Ingestion and cleaning: CSV/JSON loading, unit normalisation, quarantine.

Rules (each quarantined row records exactly one reason):
  missing_id       assay/compound/target id empty
  bad_value        value missing, non-numeric, or <= 0
  bad_unit         unit not in {nM, uM, mM}
  bad_assay_type   assay_type not in {IC50, EC50, Ki, Kd}
  out_of_range     normalised value outside [1e-3, 1e8] nM (implausible)
  duplicate        exact duplicate of an already-kept measurement
Kept rows gain ``value_nm`` and ``pic50`` (=-log10 molar potency).
"""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any, Iterable

UNIT_TO_NM = {"nm": 1.0, "um": 1e3, "µm": 1e3, "mm": 1e6}
VALID_ASSAY_TYPES = {"IC50", "EC50", "KI", "KD"}
MIN_NM, MAX_NM = 1e-3, 1e8


def load_records(path: str | Path) -> list[dict[str, Any]]:
    p = Path(path)
    if p.suffix.lower() == ".json":
        data = json.loads(p.read_text(encoding="utf-8"))
        if not isinstance(data, list):
            raise ValueError("JSON input must be a list of records")
        return data
    if p.suffix.lower() in {".csv", ".tsv"}:
        delim = "\t" if p.suffix.lower() == ".tsv" else ","
        with p.open(newline="", encoding="utf-8") as fh:
            return list(csv.DictReader(fh, delimiter=delim))
    raise ValueError(f"unsupported input format: {p.suffix}")


def write_csv(rows: Iterable[dict[str, Any]], path: str | Path, fieldnames: list[str] | None = None) -> None:
    rows = list(rows)
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        fieldnames = list(rows[0].keys()) if rows else ["assay_id"]
    with Path(path).open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def clean_records(records: Iterable[dict[str, Any]]) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, int]]:
    clean: list[dict[str, Any]] = []
    quarantine: list[dict[str, Any]] = []
    reasons: dict[str, int] = {}
    seen: set[tuple] = set()

    def reject(row: dict, reason: str) -> None:
        quarantine.append({**row, "quarantine_reason": reason})
        reasons[reason] = reasons.get(reason, 0) + 1

    for raw in records:
        row = {str(k).strip().lower(): (v.strip() if isinstance(v, str) else v) for k, v in raw.items()}
        if not row.get("assay_id") or not row.get("compound_id") or not row.get("target_id"):
            reject(row, "missing_id"); continue
        try:
            value = float(row.get("value"))
        except (TypeError, ValueError):
            reject(row, "bad_value"); continue
        if not math.isfinite(value) or value <= 0:
            reject(row, "bad_value"); continue
        unit = str(row.get("unit", "nM")).strip().lower()
        if unit not in UNIT_TO_NM:
            reject(row, "bad_unit"); continue
        assay_type = str(row.get("assay_type", "")).upper()
        if assay_type not in VALID_ASSAY_TYPES:
            reject(row, "bad_assay_type"); continue
        value_nm = value * UNIT_TO_NM[unit]
        if not (MIN_NM <= value_nm <= MAX_NM):
            reject(row, "out_of_range"); continue
        key = (row["compound_id"], row["target_id"], assay_type, round(value_nm, 6), row.get("cell_line"))
        if key in seen:
            reject(row, "duplicate"); continue
        seen.add(key)
        clean.append({
            "assay_id": str(row["assay_id"]), "compound_id": str(row["compound_id"]),
            "target_id": str(row["target_id"]), "cell_line": row.get("cell_line") or None,
            "assay_type": assay_type, "value_nm": round(value_nm, 6),
            "pic50": round(9.0 - math.log10(value_nm), 4),
            "relation": row.get("relation") or "=",
            "source": row.get("source") or "unknown",
            "is_active": (9.0 - math.log10(value_nm)) >= 6.0,
        })
    summary = {"input": len(clean) + len(quarantine), "clean": len(clean),
               "quarantined": len(quarantine), "reasons": reasons}
    return clean, quarantine, summary  # type: ignore[return-value]


def cleaning_accuracy(summary_input: int, quarantined: int, injected_defects: int) -> float:
    """Fraction of injected defects caught, assuming defects == expected quarantine.

    Exact on synthetic data: every injected defect maps to one quarantine rule
    and clean synthetic rows never trip a rule, so accuracy is 1.0 when
    quarantined == injected_defects. Returns a conservative ratio otherwise.
    """
    if injected_defects == 0:
        return 1.0 if quarantined == 0 else max(0.0, 1.0 - quarantined / max(1, summary_input))
    caught = min(quarantined, injected_defects)
    false_quarantine = max(0, quarantined - injected_defects)
    return max(0.0, (caught - false_quarantine) / injected_defects)
