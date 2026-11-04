"""Step 1 - build the dataset.

Usage (run inside the results folder):  python 01_build_dataset.py <folder with *_cells.csv>
Output: dataset.pkl = (cells, divisions, pairs), read by all later steps.
"""
import os, sys, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pandas as pd
from dataset import build, MAX_GEN

cells, div, pairs = build(sys.argv[1])
pd.to_pickle((cells, div, pairs), "dataset.pkl")
print(f"curated generations {MAX_GEN}: {len(cells):,} cells, {len(div):,} division x marker rows, {len(pairs):,} pairs")
