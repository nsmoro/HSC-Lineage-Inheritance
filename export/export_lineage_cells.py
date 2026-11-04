"""Export one compact cell table per experiment from the CURATED Matilda trees.

Usage (on the cluster):
    python export_lineage_cells.py --output ./lineage_tables                # all experiments in markers.yaml
    python export_lineage_cells.py --output ./lineage_tables --experiments 191125FS13

Reads the curated .mat files (read-only) and writes per experiment:
    <exp>_cells.csv   one row per tracked cell: lineage (Colony, TrackNumber, Generation, mother),
                      birth/last frame, fate, cell-cycle length, and for every marker its value at
                      birth (<Marker>_birth = first measurement within birth_window frames) and mean
    <exp>_qc.csv      tracking issues found while cleaning
These tables are small (a few MB) and are the input of the lineage-inheritance analysis.
"""
import argparse
import glob
import re
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from io_readers import read_position_mat, read_tat_xml
from lineage import KEY, build_cells, clean_tracked_rows

HERE = Path(__file__).resolve().parent
BASE = ["Position", "TimePoint", "Colony", "TrackNumber", "Area", "Division", "Apoptosis", "Lost"]


def summarise(df, chmap, window, issues):
    """df: rows of one file (may contain several positions). chmap: {channel: marker}."""
    for c in ["Division", "Apoptosis", "Lost"]:
        df[c] = df[c].fillna(0) if c in df else 0
    if "Area" not in df:
        df["Area"] = np.nan
    chans = []
    for ch, marker in chmap.items():
        col = f"{PREFIX}{int(ch):02d}"
        if col in df:
            df[f"Ch{int(ch):02d}"] = df.pop(col)
            chans.append(int(ch))
        else:
            issues.append({"Position": -1, "Colony": 0, "TrackNumber": 0, "Issue": "missing_channel",
                           "Detail": f"{col} ({marker}) not in file"})
    out = []
    for pos, d in df.groupby("Position"):
        rows, iss = clean_tracked_rows(d)
        issues += iss
        if rows.empty:
            continue
        cells = build_cells(rows, chans, int(d.TimePoint.max()), issues)
        # marker value at birth: first measurement within `window` frames after the cell's first frame
        r = rows.merge(cells[KEY + ["FirstTP"]], on=KEY)
        early = r[r.TimePoint < r.FirstTP + window].sort_values("TimePoint")
        for ch in chans:
            m = chmap[ch] if ch in chmap else chmap[str(ch)]
            b = early.dropna(subset=[f"Ch{ch:02d}"]).groupby(KEY)[f"Ch{ch:02d}"].first().rename(f"{m}_birth")
            cells = cells.merge(b, left_on=KEY, right_index=True, how="left")
            cells = cells.rename(columns={f"Ch{ch:02d}_Mean": f"{m}_mean"})
            cells = cells.drop(columns=[f"Ch{ch:02d}_First", f"Ch{ch:02d}_Last"])
        out.append(cells)
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()


def main():
    global PREFIX
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--output", required=True)
    ap.add_argument("--markers", default=HERE.parent / "config" / "markers.yaml")
    ap.add_argument("--experiments", nargs="*")
    a = ap.parse_args()
    cfg = yaml.safe_load(open(a.markers))
    PREFIX = cfg["intensity_prefix"]
    out = Path(a.output)
    out.mkdir(parents=True, exist_ok=True)
    roots = {"human_root": cfg.get("human_root", ""), "mouse_root": cfg.get("mouse_root", "")}

    for exp, e in cfg["experiments"].items():
        if a.experiments and exp not in a.experiments:
            continue
        f_out = out / f"{exp}_cells.csv"
        if f_out.exists():
            print(f"{exp}: {f_out} exists, skipped")
            continue
        folder = Path(e["folder"].format(**roots))
        issues, tables = [], []
        if e["species"] == "human":
            for cond, c in e["conditions"].items():
                f = folder / f"{exp}_Condition{int(cond):03d}.mat"
                if not f.exists():
                    print(f"  {exp} condition {cond}: file missing ({f.name})")
                    continue
                chmap = {int(k): v for k, v in c["markers"].items()}
                cols = BASE + [f"{PREFIX}{ch:02d}" for ch in chmap]
                t = summarise(read_position_mat(f, cols), chmap, e["birth_window"], issues)
                t.insert(0, "Condition", c["name"]); t.insert(1, "ConditionID", int(cond))
                tables.append(t)
                print(f"  {exp} condition {cond} ({c['name']}): {len(t)} cells, "
                      f"{(t.Fate == 'divided').sum()} divisions")
        else:
            _, pos_tab, _ = read_tat_xml(e["xml"].format(**roots))
            cond_of = dict(zip(pos_tab.Position, pos_tab.Condition))
            for f in sorted(folder.glob(f"{exp}_p*.mat")):
                pos = int(re.search(r"_p(\d+)", f.stem)[1])
                cname = cond_of.get(pos, "")
                panel = next((p for key, p in e["panels"].items() if key in cname), None)
                if panel is None:
                    continue
                chmap = {int(k): v for k, v in panel.items()}
                cols = BASE + [f"{PREFIX}{ch:02d}" for ch in chmap]
                df = read_position_mat(f, cols)
                if "Colony" not in df or "TrackNumber" not in df or df.empty:
                    continue
                if "Position" not in df:
                    df["Position"] = pos
                t = summarise(df, chmap, e["birth_window"], issues)
                if t.empty:
                    continue
                t.insert(0, "Condition", cname); t.insert(1, "ConditionID", 0)
                tables.append(t)
            print(f"  {exp}: {sum(len(t) for t in tables)} cells from {len(tables)} positions")
        if not tables:
            print(f"{exp}: nothing exported")
            continue
        cells = pd.concat(tables, ignore_index=True)
        cells.insert(0, "Species", e["species"]); cells.insert(0, "Experiment", exp)
        cells.round(4).to_csv(f_out, index=False)
        pd.DataFrame(issues).to_csv(out / f"{exp}_qc.csv", index=False)
        print(f"{exp}: {len(cells)} cells -> {f_out}")


if __name__ == "__main__":
    main()
