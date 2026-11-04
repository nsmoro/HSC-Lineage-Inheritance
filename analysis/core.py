"""Core data structures: divisions and mother/daughter pairs.

Every row of a cell table (<experiment>_cells.csv, written by export/export_lineage_cells.py) is one tracked cell.
TrackNumber follows the Matilda convention: the founder is 1 and the daughters of cell n are 2n and 2n+1,
so Generation = floor(log2(TrackNumber)).

This module turns cells into
  * divisions  - one row per division and marker, with four asymmetry measures between the two sisters
  * pairs      - each division joined with the division that produced its mother (same marker)

Asymmetry measures (a = sister with more marker, b = sister with less, both at their first value after birth):
  log2ratio    log2(a / b)                         -> "ratio"
  pctdiff      100 * (a - b) / mean(a, b)          -> "% difference" (same ranking as the ratio)
  absdiff      a - b                               -> "absolute difference" (depends on marker level)
  absdiff_adj  rank(absdiff) after removing its linear dependence on rank(a + b) within each stratum
For every measure we also store its percentile within experiment x condition x marker x generation
(<measure>_pct) and the corresponding thirds (<measure>_class: sym / int / asym; chance = 1/3).
"""
import glob
import numpy as np
import pandas as pd

METRICS = {"log2ratio": "Ratio (log2)", "pctdiff": "% difference", "absdiff": "Absolute difference",
           "absdiff_adj": "Abs. difference,\nadjusted for level"}
KEY = ["Experiment", "Condition", "Position", "Colony"]

def load(folder):
    """Read and concatenate all <experiment>_cells.csv tables in a folder."""
    return pd.concat([pd.read_csv(f) for f in sorted(glob.glob(f"{folder}/*_cells.csv"))], ignore_index=True)

def markers_of(cells):
    """Marker names present in a cell table (columns <marker>_birth)."""
    return [c[:-6] for c in cells.columns if c.endswith("_birth")]

def division_table(cells):
    """One row per division per marker: both daughters' birth values and the three asymmetry scores."""
    c = cells.copy()
    c["Mother"] = c.TrackNumber // 2
    rows = []
    mk = markers_of(cells)
    d2 = c[c.TrackNumber % 2 == 0].set_index(KEY + ["Mother"])
    d3 = c[(c.TrackNumber % 2 == 1) & (c.TrackNumber > 1)].set_index(KEY + ["Mother"])
    both = d2.join(d3, lsuffix="_2", rsuffix="_3", how="inner")
    base = both.reset_index()[KEY + ["Mother", "Species_2", "Generation_2", "FirstTP_2", "Fate_2", "Fate_3"]]
    base = base.rename(columns={"Species_2": "Species", "Generation_2": "DaughterGen", "FirstTP_2": "BirthTP"})
    for m in mk:
        a, b = both[f"{m}_birth_2"].values, both[f"{m}_birth_3"].values
        ok = ~np.isnan(a) & ~np.isnan(b)
        hi, lo = np.maximum(a, b), np.minimum(a, b)
        with np.errstate(divide="ignore", invalid="ignore"):
            t = base.copy()
            t["Marker"] = m
            t["high_daughter"] = np.where(a >= b, 2, 3)          # 2 = even TrackNumber (2n), 3 = odd (2n+1)
            t["v_high"], t["v_low"] = hi, lo
            t["log2ratio"] = np.where((lo > 0) & ok, np.log2(hi / lo), np.nan)
            t["absdiff"] = np.where(ok, hi - lo, np.nan)
            t["pctdiff"] = np.where(ok & (hi + lo > 0), 100 * (hi - lo) / ((hi + lo) / 2), np.nan)
        rows.append(t[ok])
    d = pd.concat(rows, ignore_index=True)
    # percentile of each score within its stratum (experiment x condition x marker x generation):
    # removes differences in level between conditions and generations, keeps the ranking
    strata = ["Experiment", "Condition", "Marker", "DaughterGen"]
    # absolute difference adjusted for the amount of marker (sum of both daughters):
    # residual of rank(absdiff) after a linear fit on rank(total), within each stratum
    d["total"] = d.v_high + d.v_low
    def _adj(g):
        y = g.absdiff.rank(pct=True); x = g.total.rank(pct=True)
        if len(g) < 5 or x.nunique() < 2:
            return pd.Series(np.nan, index=g.index)
        b = np.polyfit(x, y, 1)
        return y - np.polyval(b, x)
    d["absdiff_adj"] = d.groupby(strata, group_keys=False).apply(_adj)
    for s in METRICS:
        d[f"{s}_pct"] = d.groupby(strata)[s].rank(pct=True)
        d[f"{s}_class"] = pd.cut(d[f"{s}_pct"], [0, 1/3, 2/3, 1.0001], labels=["sym", "int", "asym"])
    return d

def parent_child(div):
    """Pair each division with the division that produced its mother (same marker)."""
    p = div.copy()
    child = p.rename(columns=lambda c: c if c in KEY + ["Marker"] else "c_" + c)
    parent = p.rename(columns=lambda c: c if c in KEY + ["Marker"] else "p_" + c)
    child["p_Mother"] = child["c_Mother"] // 2
    pairs = child.merge(parent, on=KEY + ["Marker", "p_Mother"])
    # did the dividing mother receive the higher amount in the parent division?
    pairs["mother_was_high"] = (pairs.c_Mother % 2 == 0) == (pairs.p_high_daughter == 2)
    return pairs
