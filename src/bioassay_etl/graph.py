"""Transformation to a knowledge-graph structure (stdlib, no graph DB needed).

Node types : Compound, Target, Assay, CellLine
Edge types : TESTED_IN   Compound -> Assay
             MEASURES    Assay -> Target
             IN_CELL     Assay -> CellLine        (when a cell line is present)
             ACTIVE_AGAINST Compound -> Target    (derived, mean pIC50 >= 6.0)
The graph is exported as nodes.csv / edges.csv / graph.json so it can be
loaded into Neo4j, NetworkX, or a dataframe tool without this package.
"""
from __future__ import annotations

import json
from collections import defaultdict, deque
from pathlib import Path
from typing import Any, Iterable

ACTIVE_PIC50 = 6.0


class KnowledgeGraph:
    def __init__(self) -> None:
        self.nodes: dict[str, dict[str, Any]] = {}   # id -> {id, type, **props}
        self.edges: list[dict[str, Any]] = []
        self._adj: dict[str, list[str]] = defaultdict(list)
        self._edge_keys: set[tuple[str, str, str]] = set()

    def add_node(self, node_id: str, node_type: str, **props: Any) -> None:
        if node_id not in self.nodes:
            self.nodes[node_id] = {"id": node_id, "type": node_type, **props}

    def add_edge(self, src: str, dst: str, edge_type: str, **props: Any) -> None:
        key = (src, dst, edge_type)
        if key in self._edge_keys:
            return
        self._edge_keys.add(key)
        self.edges.append({"src": src, "dst": dst, "type": edge_type, **props})
        self._adj[src].append(dst)
        self._adj[dst].append(src)  # undirected traversal for neighbourhood queries

    # -- queries ---------------------------------------------------------
    def neighbors(self, node_id: str, depth: int = 1) -> set[str]:
        seen, frontier = {node_id}, {node_id}
        for _ in range(depth):
            nxt = {m for f in frontier for m in self._adj.get(f, [])} - seen
            seen |= nxt
            frontier = nxt
        return seen - {node_id}

    def degree(self, node_id: str) -> int:
        return len(self._adj.get(node_id, []))

    def top_hubs(self, k: int = 10, node_type: str | None = None) -> list[tuple[str, int]]:
        items = [(nid, self.degree(nid)) for nid, n in self.nodes.items()
                 if node_type is None or n["type"] == node_type]
        return sorted(items, key=lambda x: -x[1])[:k]

    def stats(self) -> dict[str, Any]:
        by_type: dict[str, int] = defaultdict(int)
        for n in self.nodes.values():
            by_type[n["type"]] += 1
        by_edge: dict[str, int] = defaultdict(int)
        for e in self.edges:
            by_edge[e["type"]] += 1
        return {"nodes": len(self.nodes), "edges": len(self.edges),
                "nodes_by_type": dict(by_type), "edges_by_type": dict(by_edge)}

    # -- export ----------------------------------------------------------
    def to_dict(self) -> dict[str, Any]:
        return {"nodes": list(self.nodes.values()), "edges": self.edges}

    def export(self, out_dir: str | Path) -> dict[str, str]:
        from .ingest import write_csv
        out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)
        write_csv(self.nodes.values(), out / "nodes.csv", ["id", "type"])
        write_csv(self.edges, out / "edges.csv", ["src", "dst", "type", "pic50", "assay_type"])
        (out / "graph.json").write_text(json.dumps(self.to_dict(), indent=1), encoding="utf-8")
        return {"nodes": str(out / "nodes.csv"), "edges": str(out / "edges.csv"),
                "graph_json": str(out / "graph.json")}


def build_graph(records: Iterable[dict[str, Any]]) -> KnowledgeGraph:
    g = KnowledgeGraph()
    potency: dict[tuple[str, str], list[float]] = defaultdict(list)
    for r in records:
        cid, tid, aid = f"compound:{r['compound_id']}", f"target:{r['target_id']}", f"assay:{r['assay_id']}"
        g.add_node(cid, "Compound")
        g.add_node(tid, "Target")
        g.add_node(aid, "Assay", assay_type=r["assay_type"], pic50=r["pic50"], value_nm=r["value_nm"])
        g.add_edge(cid, aid, "TESTED_IN")
        g.add_edge(aid, tid, "MEASURES", pic50=r["pic50"], assay_type=r["assay_type"])
        if r.get("cell_line"):
            cell = f"cell:{r['cell_line']}"
            g.add_node(cell, "CellLine")
            g.add_edge(aid, cell, "IN_CELL")
        potency[(cid, tid)].append(float(r["pic50"]))
    for (cid, tid), vals in potency.items():
        mean_pic50 = sum(vals) / len(vals)
        if mean_pic50 >= ACTIVE_PIC50:
            g.add_edge(cid, tid, "ACTIVE_AGAINST", pic50=round(mean_pic50, 4))
    return g
