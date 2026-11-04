"""Restrict the cell tables to the curated generations and build the analysis dataset.

Curated generations (see README): mouse cells of generations 0-3, human all curated cells (effectively
generations 0-2; a few generation-3 cells exist in one experiment). The mouse limit can be changed with the
environment variable MAX_GEN_MOUSE (used for the generation-4 robustness check).
"""
import os
import numpy as np
import pandas as pd
from core import load, division_table, parent_child, KEY  # noqa: F401  (KEY re-exported for the scripts)

MAX_GEN = {"mouse": int(os.environ.get("MAX_GEN_MOUSE", 3)), "human": 3}


def build(tables_folder):
    """Return (cells, divisions, mother/daughter pairs) for the curated generations."""
    cells = load(tables_folder)
    cells = cells[cells.Generation <= cells.Species.map(MAX_GEN)].copy()
    div = division_table(cells)
    pairs = parent_child(div)
    return cells, div, pairs
