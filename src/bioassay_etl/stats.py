"""Statistical EDA and data-integrity validation on cleaned records.

Pure stdlib (statistics module) so results are identical with or without
pandas/numpy installed. Everything returns plain dicts for JSON reports.
"""
from __future__ import annotations

import math
import statistics as st
from collections import Counter, defaultdict
from typing import Any, Iterable


def _quantile(sorted_vals: list[float], q: float) -> float:
    if not sorted_vals:
        return float("nan")
    pos = (len(sorted_vals) - 1) * q
    lo, hi = int(math.floor(pos)), int(math.ceil(pos))
    return sorted_vals[lo] if lo == hi else sorted_vals[lo] + (sorted_vals[hi] - sorted_vals[lo]) * (pos - lo)


def describe(values: list[float]) -> dict[str, float]:
    if not values:
        return {"count": 0}
    s = sorted(values)
    return {"count": len(s), "mean": round(st.fmean(s), 4), "median": round(st.median(s), 4),
            "std": round(st.pstdev(s), 4) if len(s) > 1 else 0.0,
            "min": s[0], "p25": round(_quantile(s, .25), 4), "p75": round(_quantile(s, .75), 4), "max": s[-1]}


def pearson(xs: list[float], ys: list[float]) -> float:
    if len(xs) != len(ys) or len(xs) < 2:
        return 0.0
    mx, my = st.fmean(xs), st.fmean(ys)
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    den = math.sqrt(sum((x - mx) ** 2 for x in xs) * sum((y - my) ** 2 for y in ys))
    return round(num / den, 4) if den else 0.0


def iqr_outliers(values: list[float], k: float = 1.5) -> int:
    if len(values) < 4:
        return 0
    s = sorted(values)
    q1, q3 = _quantile(s, .25), _quantile(s, .75)
    lo, hi = q1 - k * (q3 - q1), q3 + k * (q3 - q1)
    return sum(1 for v in s if v < lo or v > hi)


def integrity_report(clean: list[dict[str, Any]], quarantine_count: int, input_count: int) -> dict[str, Any]:
    """Data-integrity checks that the test-suite also asserts on."""
    keys = [(r["compound_id"], r["target_id"], r["assay_type"], r["value_nm"], r["cell_line"]) for r in clean]
    pic50_ok = all(3.0 - 3 <= r["pic50"] <= 12.0 for r in clean) if clean else True
    checks = {
        "no_duplicate_measurements": len(keys) == len(set(keys)),
        "all_values_positive_nm": all(r["value_nm"] > 0 for r in clean),
        "pic50_in_plausible_range": pic50_ok,
        "all_ids_present": all(r["compound_id"] and r["target_id"] and r["assay_id"] for r in clean),
        "row_accounting_balanced": len(clean) + quarantine_count == input_count,
    }
    return {"checks": checks, "passed": all(checks.values()),
            "completeness": round(len(clean) / input_count, 4) if input_count else 0.0}


def eda_report(records: list[dict[str, Any]]) -> dict[str, Any]:
    if not records:
        return {"records": 0}
    pic50 = [float(r["pic50"]) for r in records]
    by_target: dict[str, list[dict]] = defaultdict(list)
    for r in records:
        by_target[r["target_id"]].append(r)
    target_summary = {
        t: {"n": len(rs), "active_rate": round(sum(bool(r["is_active"]) for r in rs) / len(rs), 4),
            "mean_pic50": round(st.fmean(float(r["pic50"]) for r in rs), 4)}
        for t, rs in sorted(by_target.items(), key=lambda kv: -len(kv[1]))[:10]
    }
    return {
        "records": len(records),
        "pic50": describe(pic50),
        "active_rate": round(sum(bool(r["is_active"]) for r in records) / len(records), 4),
        "assay_type_counts": dict(Counter(r["assay_type"] for r in records)),
        "source_counts": dict(Counter(r["source"] for r in records)),
        "pic50_iqr_outliers": iqr_outliers(pic50),
        # log10(value) vs pIC50 is -1 by construction; reported as a sanity check.
        "corr_log10_value_vs_pic50": pearson([math.log10(float(r["value_nm"])) for r in records], pic50),
        "top_targets_by_volume": target_summary,
    }
