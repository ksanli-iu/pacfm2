"""PACFM 2 — specificity-aware pathway profiling of KEGG orthologue (KO) abundance tables.

Most KOs take part in several KEGG pathways, so summing every KO of a pathway ("naive" profiling) attributes the
abundance of widely shared enzymes to every pathway that contains them, including pathways that are absent from the
sample (e.g. human disease pathways in bacterial communities). This module computes pathway abundances under several
counting schemes and, for every pathway and sample, the evidence that comes from enzymes specific to that pathway.

Schemes
  naive        sum of all KOs of the pathway
  cutoff       PACFM 2016: only KOs associated with <= c pathways (c = 1: pathway-unique KOs only)
  fractional   each KO's abundance divided equally among its pathways (abundance is conserved)
  em           evidence-proportional allocation: a shared KO's abundance is divided among its pathways in proportion
               to their estimated per-KO abundance, re-estimated iteratively (EM-like; abundance is conserved)
  fraccov      coverage-weighted fractional: a shared KO is divided among its pathways in proportion to each pathway's
               KO coverage in the sample (detected KOs / pathway KOs), raised to the power gamma (default 2)
  emcov        em with each pathway's share additionally weighted by its KO coverage ** gamma (regularised em: a nearly
               complete pathway keeps its share even when all its KOs are shared)

Global and overview maps (KEGG 011xx/012xx) are excluded, since they are unions of other pathways.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

GLOBAL_MAPS = re.compile(r"^01[12]\d\d$")


# ------------------------------------------------------------------------------------------------ KEGG reference
@dataclass
class Kegg:
    ko_pathway: pd.DataFrame                 # columns ko, pathway
    hierarchy: pd.DataFrame                  # pathway -> level_a, level_b, name
    ko_module: pd.DataFrame | None = None
    module_def: pd.DataFrame | None = None   # module -> name, class, definition
    module_pathway: pd.DataFrame | None = None
    n_assoc: pd.Series = field(init=False)   # KO -> number of pathways

    def __post_init__(self):
        kp = self.ko_pathway
        kp = kp[~kp.pathway.str.match(GLOBAL_MAPS)]
        kp = kp[kp.pathway.isin(self.hierarchy.pathway)].drop_duplicates()
        self.ko_pathway = kp.reset_index(drop=True)
        self.n_assoc = kp.groupby("ko").pathway.nunique()

    @classmethod
    def load(cls, d: str | Path) -> "Kegg":
        d = Path(d)
        rd = lambda f: pd.read_csv(d / f, sep="\t", dtype=str) if (d / f).exists() else None
        return cls(rd("ko_pathway.tsv"), rd("pathway_hierarchy.tsv"), rd("ko_module.tsv"), rd("module_definition.tsv"),
                   rd("module_pathway.tsv"))

    def pathway_kos(self) -> dict[str, set[str]]:
        return self.ko_pathway.groupby("pathway").ko.apply(set).to_dict()

    def specificity(self) -> pd.DataFrame:
        """per pathway: KOs, KOs unique to the pathway, median number of pathways per KO"""
        kp = self.ko_pathway.assign(n=lambda x: x.ko.map(self.n_assoc))
        s = kp.groupby("pathway").agg(kos=("ko", "size"), unique_kos=("n", lambda v: int((v == 1).sum())),
                                      median_assoc=("n", "median"))
        s["pct_unique"] = 100 * s.unique_kos / s.kos
        return s.join(self.hierarchy.set_index("pathway"))


# ------------------------------------------------------------------------------------------------ profiling
def _incidence(kegg: Kegg, kos: pd.Index, max_assoc: int | None = None):
    kp = kegg.ko_pathway[kegg.ko_pathway.ko.isin(kos)]
    if max_assoc is not None:
        kp = kp[kp.ko.map(kegg.n_assoc) <= max_assoc]
    M = pd.crosstab(kp.ko, kp.pathway).clip(upper=1).astype(float)
    return M.reindex(index=kos, fill_value=0.0)


def pathway_abundance(ko_table: pd.DataFrame, kegg: Kegg, scheme: str = "naive", cutoff: int = 1,
                      max_iter: int = 200, tol: float = 1e-8, gamma: float = 2.0) -> pd.DataFrame:
    """ko_table: KO x sample abundances (index = KO ids). Returns pathway x sample abundances."""
    X = ko_table.groupby(level=0).sum()
    X = X[X.index.isin(kegg.n_assoc.index)]
    if scheme == "naive":
        M = _incidence(kegg, X.index)
        return M.T @ X
    if scheme == "cutoff":
        M = _incidence(kegg, X.index, max_assoc=cutoff)
        return M.T @ X
    M = _incidence(kegg, X.index)
    if scheme == "fractional":
        W = M.div(M.sum(axis=1), axis=0)
        return W.T @ X
    if scheme in ("em", "emcov", "fraccov"):
        size = pd.Series(kegg.ko_pathway.groupby("pathway").ko.nunique()).reindex(M.columns).values
        out = {}
        Mv = M.values
        for s in X.columns:
            a = X[s].values
            cov = (Mv.T @ (a > 0).astype(float)) / size     # KO coverage of each pathway in this sample
            prior = cov ** gamma if scheme in ("emcov", "fraccov") else np.ones_like(cov)
            if scheme == "fraccov":
                share = Mv * prior
                denom = share.sum(axis=1, keepdims=True)
                share = np.divide(share, denom, out=np.zeros_like(share), where=denom > 0)
                out[s] = share.T @ a
                continue
            theta = (Mv.T @ a) / size                      # start: naive per-KO abundance
            for _ in range(max_iter):
                share = Mv * theta * prior                  # KO x pathway weights
                denom = share.sum(axis=1, keepdims=True)
                share = np.divide(share, denom, out=np.zeros_like(share), where=denom > 0)
                total = share.T @ a
                new = total / size
                if np.max(np.abs(new - theta)) <= tol * max(1.0, np.max(theta)):
                    theta = new; break
                theta = new
            out[s] = total
        return pd.DataFrame(out, index=M.columns)
    raise ValueError(f"unknown scheme {scheme!r}")


def pathway_evidence(ko_table: pd.DataFrame, kegg: Kegg, detect: float = 0.0) -> dict[str, pd.DataFrame]:
    """per pathway x sample: coverage (detected KOs / pathway KOs), unique_coverage (detected unique KOs / unique KOs),
    unique_share (share of naive abundance carried by pathway-unique KOs)"""
    X = ko_table.groupby(level=0).sum()
    X = X[X.index.isin(kegg.n_assoc.index)]
    M = _incidence(kegg, X.index)
    U = M.mul((kegg.n_assoc.reindex(M.index) == 1).astype(float), axis=0)
    size = kegg.ko_pathway.groupby("pathway").ko.nunique().reindex(M.columns)
    usize = kegg.ko_pathway.assign(n=lambda x: x.ko.map(kegg.n_assoc)).query("n == 1").groupby("pathway").ko.nunique().reindex(M.columns)
    det = (X > detect).astype(float)
    cov = (M.T @ det).div(size, axis=0)
    ucov = (U.T @ det).div(usize, axis=0)
    naive = M.T @ X
    ushare = (U.T @ X).div(naive.replace(0, np.nan))
    return {"coverage": cov, "unique_coverage": ucov, "unique_share": ushare}


# ------------------------------------------------------------------------------------------------ KEGG modules
_TOK = re.compile(r"K\d{5}|M\d{5}|--|[(),+\- ]")


class ModuleDefinition:
    """KEGG module DEFINITION: space = consecutive steps (AND), ',' = alternatives (OR), '+' = complex subunits (AND),
    '-' = optional subunit (ignored), '--' = no KO assigned to the step (treated as satisfied)."""

    def __init__(self, definition: str, resolve=None):
        self.toks = [t for t in _TOK.findall(definition)]
        self.resolve = resolve                     # module id -> bool (nested modules)

    def steps(self) -> list[list]:
        """top-level steps as token lists"""
        out, cur, depth = [], [], 0
        for t in self.toks:
            if t == "(":
                depth += 1
            elif t == ")":
                depth -= 1
            if t == " " and depth == 0:
                if cur:
                    out.append(cur); cur = []
                continue
            cur.append(t)
        if cur:
            out.append(cur)
        return out

    def _eval(self, toks, present):
        pos = 0

        def alt():                       # alternatives separated by ','
            nonlocal pos
            v = seq()
            while pos < len(toks) and toks[pos] == ",":
                pos += 1; v = seq() or v
            return v

        def seq():                       # AND of space-separated items
            nonlocal pos
            v = comp()
            while pos < len(toks) and toks[pos] == " ":
                pos += 1; v = comp() and v
            return v

        def comp():                      # complex: a+b (required), a-b (b optional)
            nonlocal pos
            v = atom()
            while pos < len(toks) and toks[pos] in "+-" and toks[pos] != "--":
                op = toks[pos]; pos += 1
                w = atom()
                if op == "+":
                    v = v and w
            return v

        def atom():
            nonlocal pos
            t = toks[pos]; pos += 1
            if t == "(":
                v = alt(); pos += 1      # ')'
                return v
            if t == "--":
                return True
            if t.startswith("M"):
                return bool(self.resolve(t)) if self.resolve else False
            return t in present
        return alt()

    def completeness(self, present: set[str]) -> float:
        st = self.steps()
        return sum(bool(self._eval(s, present)) for s in st) / len(st) if st else 0.0


def module_completeness(kegg: Kegg, present: set[str]) -> pd.Series:
    defs = kegg.module_def.set_index("module").definition
    cache = {}

    def resolve(m):
        if m not in cache:
            cache[m] = False             # guard against cycles
            cache[m] = ModuleDefinition(defs.get(m, ""), resolve).completeness(present) == 1.0
        return cache[m]
    return pd.Series({m: ModuleDefinition(d, resolve).completeness(present) for m, d in defs.items()})
