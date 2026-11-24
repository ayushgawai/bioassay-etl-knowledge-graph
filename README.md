# bioassay-etl-knowledge-graph

End-to-end ETL pipeline and graph model for pharmaceutical bioassay data,
featuring Python (with optional PySpark and pandas paths). It ingests raw
assay measurements, cleans and quarantines defective rows with explicit
reasons, normalises potency units, transforms the result into a
compound–target knowledge graph, and validates it statistically — with
data-integrity checks as tests and a benchmark against *known* injected
defects, so accuracy is measured, not asserted.

## Design

```
                 +----------------------+
  raw CSV/JSON   |  ingest.py           |
  (or synthetic  |  load -> clean ->    |-- quarantine.csv (row + reason)
   generator) -->|  unit normalisation  |
                 |  dedupe, range check |
                 +----------+-----------+
                            | clean records (value_nm, pIC50, is_active)
                            v
                 +----------------------+      +----------------------+
                 |  graph.py            |      |  stats.py            |
                 |  Compound / Target / |      |  describe, Pearson,  |
                 |  Assay / CellLine    |      |  IQR outliers,       |
                 |  nodes + typed edges |      |  integrity checks    |
                 +----------+-----------+      +----------+-----------+
                            |                             |
                            v                             v
              nodes.csv / edges.csv / graph.json     report.json (cleaning,
                                                     integrity, EDA, graph)
  optional: spark_job.py runs the (compound, target) aggregation on PySpark
  when installed; local_activity_counts() is the tested stdlib equivalent.
```

| Module | Responsibility | Key choices |
|---|---|---|
| `src/bioassay_etl/synthetic.py` | Deterministic synthetic assays + known-defect injection | Defects injected in two phases (mutations, then duplicates of valid rows) so ground truth is exact |
| `src/bioassay_etl/ingest.py` | Load CSV/TSV/JSON, clean, quarantine | Units nM/µM/mM → `value_nm`; `pIC50 = 9 − log10(value_nm)`; one explicit reason per rejected row; exact-duplicate removal |
| `src/bioassay_etl/graph.py` | Knowledge-graph transform | Node types Compound/Target/Assay/CellLine; edges `TESTED_IN`, `MEASURES`, `IN_CELL`, derived `ACTIVE_AGAINST` (mean pIC50 ≥ 6.0); exports CSV + JSON for Neo4j/NetworkX |
| `src/bioassay_etl/stats.py` | Statistical EDA + integrity | Stdlib-only describe/quantiles/Pearson/IQR outliers; integrity: no duplicates, positive values, plausible pIC50, balanced row accounting |
| `src/bioassay_etl/spark_job.py` | Optional distributed aggregation | PySpark `(compound, target)` counts when `.[spark]` is installed; raises rather than faking it otherwise |
| `src/bioassay_etl/pipeline.py` | Orchestration | ingest → clean → graph → EDA → `clean.csv`, `quarantine.csv`, graph export, `report.json` |
| `src/bioassay_etl/cli.py` | CLI runner | `generate` / `run` / `demo` / `bench` subcommands |

## Run

```bash
pip install -e .            # core has zero dependencies; .[spark] / .[pandas] are optional
python -m bioassay_etl.cli demo --rows 5000 --out out/
python -m bioassay_etl.cli generate --rows 10000 --out data/raw.csv
python -m bioassay_etl.cli run --input data/raw.csv --out out/
# outputs: out/clean.csv, out/quarantine.csv, out/graph/{nodes,edges}.csv,
#          out/graph/graph.json, out/report.json
```

## Tests

```bash
pip install pytest
python -m pytest -q
```

Sanity run (actual output):

```
............                                                             [100%]
12 passed in 0.10s
```

Tests cover every quarantine rule, unit conversion and pIC50 maths, graph
node/edge shape and the derived `ACTIVE_AGAINST` edge, integrity-report
pass/fail, the stdlib activity aggregation, and CLI end-to-end runs.

## Benchmark

Accuracy is scored against defects the generator injected on purpose
(duplicates, bad units, missing ids, non-positive values; 0.2% defect rate):
a defect is "caught" iff its row lands in quarantine, and clean synthetic
rows never trip a rule, so the score is exact *for this synthetic data*.

```bash
python scripts/benchmark.py --rows 200000     # scaled-down sanity run
python scripts/benchmark.py --rows 1000000    # full run
```

Actual outputs from a local run (Python 3.12, Windows, single process,
pure stdlib — no Spark, no pandas):

```
$ python scripts/benchmark.py --rows 200000
rows_input=200100 clean=199700 quarantined=400 injected_defects=400
cleaning_accuracy=100.0000% nodes=239501 edges=725378
generate=0.49s clean=0.50s graph=1.37s total=2.36s
throughput=84,838 rows/s

$ python scripts/benchmark.py --rows 1000000
rows_input=1000500 clean=998500 quarantined=2000 injected_defects=2000
cleaning_accuracy=100.0000% nodes=1197172 edges=3625106
generate=2.65s clean=3.09s graph=10.52s total=16.27s
throughput=61,500 rows/s
```

So the 1M+ record scale is real and reproducible here, and cleaning caught
100% of the 2,000 injected defects in that run (≥ 99.9% target met **on
synthetic data with these four defect types** — see Limitations).

## Limitations (honest)

- **Accuracy claim is synthetic-only.** 100% / ≥99.9% means "caught the
  defects this generator injects". Real assay feeds contain messier faults
  (unit ambiguity, swapped columns, censored `<`/`>` values that bias
  pIC50, plate effects) that this rule set would not score 100% on.
- **The benchmark is single-process pure Python.** ~61k rows/s at 1M rows
  is a laptop figure, not a cluster figure. The PySpark module exists for
  the distributed aggregation but was **not** benchmarked here (PySpark is
  not installed in the dev environment), and it is not required to run
  anything else.
- **Graph is in-memory and assay-centric**: every assay row becomes a node
  (1.2M nodes at 1M rows). Fine for export/analysis; not a replacement for
  a graph database at much larger scale. Export targets Neo4j/NetworkX.
- **Cleaning rules are deliberately simple**: fixed plausible-range window
  (10⁻³–10⁸ nM), exact-duplicate matching only (no fuzzy/near-duplicate
  detection), `relation` (`<`, `>`) is carried through but not used to
  adjust potency statistics.
- **EDA is descriptive, not inferential**: no multiple-testing correction,
  no dose–response curve fitting, no target-family statistics.
- pIC50 is computed identically for IC50/EC50/Ki/Kd for simplicity;
  mixing assay types in one `ACTIVE_AGAINST` edge is a modelling shortcut.

## License

MIT — see `LICENSE`.
