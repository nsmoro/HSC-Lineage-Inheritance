"""Step 5a - Rebuild lineage trees from a Position table and check tracking quality.

TrackNumber convention (same as Matilda / the MATLAB sister script):
    founder cell = 1, daughters of cell n = 2n and 2n+1, generation = floor(log2(n))

Main functions
- build_cells(rows, ...)   -> one row per tracked cell (birth, division, fate, intensities)
- tree_summary(cells, ...) -> one row per colony (depth, completeness, fates)
"""
import numpy as np
import pandas as pd

KEY = ["Position", "Colony", "TrackNumber"]


def _issue(issues, pos, col, track, kind, detail):
    issues.append({"Position": pos, "Colony": col, "TrackNumber": track,
                   "Issue": kind, "Detail": detail})


def clean_tracked_rows(df, max_timepoint=None):
    """Keep tracked rows only and fix/record the tracking quirks.

    Returns (clean_rows, issues) where issues is a list of dicts for the QC report.
    """
    issues = []
    t = df[(df["Colony"] >= 1) & (df["TrackNumber"] >= 1)].copy()
    n_colony0 = int(((df["Colony"] == 0) & (df["TrackNumber"] >= 1)).sum())
    if n_colony0:
        pos = int(df["Position"].iloc[0])
        _issue(issues, pos, 0, 0, "track_without_colony",
               f"{n_colony0} rows have a TrackNumber but Colony 0 (ignored)")
    if max_timepoint:
        t = t[t["TimePoint"] <= max_timepoint]
    t = t.sort_values(KEY + ["TimePoint"])

    # 1) duplicate rows: same cell twice at the same timepoint -> keep the first
    dup = t.duplicated(KEY + ["TimePoint"], keep="first")
    for key, n in t[dup].groupby(KEY).size().items():
        _issue(issues, *key, "duplicate_rows", f"{n} duplicated timepoints removed")
    t = t[~dup]

    # 2) rows labelled with a cell after its daughters already appeared -> removed
    first = t.groupby(KEY)["TimePoint"].min().rename("FirstTP").reset_index()
    kids = first.copy()
    kids["TrackNumber"] = kids["TrackNumber"] // 2          # daughter -> its mother
    kids = kids.groupby(KEY)["FirstTP"].min().rename("DaughtersFirstTP").reset_index()
    t = t.merge(kids, on=KEY, how="left")
    stray = t["DaughtersFirstTP"].notna() & (t["TimePoint"] >= t["DaughtersFirstTP"])
    for key, g in t[stray].groupby(KEY):
        _issue(issues, *key, "rows_after_division",
               f"{len(g)} rows at TP {g.TimePoint.min()}-{g.TimePoint.max()} "
               f"after daughters appeared (removed)")
    t = t[~stray].drop(columns="DaughtersFirstTP")
    return t, issues


def build_cells(rows, channels, movie_last_tp, issues):
    """Summarise each tracked cell. `rows` must come from clean_tracked_rows."""
    g = rows.groupby(KEY)
    cells = g.agg(FirstTP=("TimePoint", "min"), LastTP=("TimePoint", "max"),
                  nRows=("TimePoint", "size"),
                  nDivisionFlags=("Division", "sum"),
                  ApoptosisFlag=("Apoptosis", "max"),
                  LostFlag=("Lost", "max")).reset_index()
    cells["Generation"] = np.floor(np.log2(cells["TrackNumber"])).astype(int)
    cells["MotherTrack"] = cells["TrackNumber"] // 2
    cells["Area_First"] = g["Area"].first().values

    # intensities: first available value, last value, mean over the cell's life
    for ch in channels:
        col = f"Ch{ch:02d}"
        cells[f"{col}_First"] = g[col].first().values   # first non-missing value
        cells[f"{col}_Last"] = g[col].last().values
        cells[f"{col}_Mean"] = g[col].mean().values

    # which daughters exist?
    present = set(zip(cells.Position, cells.Colony, cells.TrackNumber))
    has_d2 = [(p, c, 2 * n) in present for p, c, n in present_iter(cells)]
    has_d3 = [(p, c, 2 * n + 1) in present for p, c, n in present_iter(cells)]
    cells["nDaughtersTracked"] = np.array(has_d2, int) + np.array(has_d3, int)

    # fate of each cell
    fate = np.select(
        [cells.nDaughtersTracked > 0,
         cells.ApoptosisFlag > 0,
         cells.LostFlag > 0,
         cells.LastTP >= movie_last_tp,
         cells.nDivisionFlags > 0],
        ["divided", "died", "lost", "end_of_movie", "divided_daughters_not_tracked"],
        default="track_ends")
    cells["Fate"] = fate
    cells["CellCycleTPs"] = np.where((cells.Fate == "divided") & (cells.TrackNumber > 1),
                                     cells.LastTP - cells.FirstTP + 1, np.nan)

    # quality checks
    for r in cells.itertuples():
        key = (r.Position, r.Colony, r.TrackNumber)
        span = r.LastTP - r.FirstTP + 1
        if r.TrackNumber > 1 and (r.Position, r.Colony, r.MotherTrack) not in present:
            _issue(issues, *key, "mother_missing", "cell has no tracked mother")
        if r.nDaughtersTracked == 1:
            _issue(issues, *key, "one_daughter_only", "only one of two daughters tracked")
        if r.nDivisionFlags > 1:
            _issue(issues, *key, "repeated_division_flag",
                   f"Division flag set on {int(r.nDivisionFlags)} rows")
        if r.nDaughtersTracked > 0 and r.nDivisionFlags == 0:
            _issue(issues, *key, "daughters_without_division_flag",
                   "daughters exist but no Division flag")
        if r.nRows < span:
            _issue(issues, *key, "gaps", f"{span - r.nRows} of {span} timepoints missing")
    return cells


def present_iter(cells):
    return zip(cells.Position, cells.Colony, cells.TrackNumber)


def tree_summary(cells, sister_pairs=None):
    """One row per colony: how deep and how completely it was tracked."""
    out = []
    for (pos, col), c in cells.groupby(["Position", "Colony"]):
        fate = dict(zip(c.TrackNumber, c.Fate))
        daughters = set(c.TrackNumber)
        # complete generation k: every cell of generation k-1 either divided with both
        # daughters tracked, or died
        complete, parents = 0, [1] if 1 in fate else []
        while parents:
            nxt, ok = [], True
            for n in parents:
                if fate[n] == "died":
                    continue
                if 2 * n in daughters and 2 * n + 1 in daughters:
                    nxt += [2 * n, 2 * n + 1]
                else:
                    ok = False
                    break
            if not ok or not nxt:
                break
            complete += 1
            parents = nxt
        row = {"Condition": c.Condition.iloc[0], "Position": pos, "Colony": col,
               "nCells": len(c), "MaxGeneration": int(c.Generation.max()),
               "CompleteGenerations": complete,
               "nDivided": int((c.Fate == "divided").sum()),
               "nDied": int((c.Fate == "died").sum()),
               "nLost": int((c.Fate == "lost").sum()),
               "nTrackEnds": int(c.Fate.isin(["track_ends",
                                              "divided_daughters_not_tracked"]).sum()),
               "nEndOfMovie": int((c.Fate == "end_of_movie").sum())}
        if sister_pairs is not None and len(sister_pairs):
            s = sister_pairs[(sister_pairs.Position == pos) & (sister_pairs.Colony == col)]
            for cls in ["symmetric", "intermediate", "asymmetric"]:
                row[f"n_{cls}"] = int((s.Class == cls).sum())
        out.append(row)
    return pd.DataFrame(out).sort_values(["CompleteGenerations", "MaxGeneration", "nCells"],
                                         ascending=False)
