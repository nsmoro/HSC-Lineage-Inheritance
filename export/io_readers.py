"""Step 1 - Reading the input files (read-only: nothing here writes to disk).

- read_tat_xml(path)            -> experiment info, positions table, channels table
- read_position_mat(path, cols) -> one Matilda Position .mat file as a pandas table
"""
import re
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io as sio


def read_tat_xml(path):
    """Read the microscope metadata file <experiment>_TATexp.xml."""
    root = ET.parse(path).getroot()
    exp_name = Path(path).name.split("_TATexp")[0]

    positions = []
    for p in root.iter("PosInfoDimension"):
        comment = p.get("comments", "").strip()
        # comment looks like "A1 [1:1] HSC Replicate 1 FcyRBV480 ..." -> drop the well label
        condition = re.sub(r"^\S+\s+\[\d+:\d+\]\s*", "", comment).strip()
        positions.append({"Experiment": exp_name, "Position": int(p.get("index")),
                          "Well": comment.split(" ")[0] if comment else "",
                          "Condition": condition})

    channels = []
    for w in root.iter("WLInfo"):
        channels.append({"Channel": int(w.get("Name")), "Marker": w.get("Comment"),
                         "ExposureMs": float(w.get("ExpTime")),
                         "EveryNthTimepoint": int(w.get("Interval")),
                         "ImageType": w.get("ImageType")})

    mpp = root.find(".//MicrometerPerPixel")
    meta = {"Experiment": exp_name,
            "MicrometerPerPixel": float(mpp.get("value")) if mpp is not None else None}
    return meta, pd.DataFrame(positions), pd.DataFrame(channels)


def read_position_mat(path, columns=None):
    """Read one Position .mat file (one row per segmented object per timepoint).

    columns: list of variable names to load (loading only what is needed saves memory).
    Variables missing from the file are skipped.
    """
    m = sio.loadmat(path, squeeze_me=True, variable_names=columns)
    data = {k: np.atleast_1d(v) for k, v in m.items() if not k.startswith("__")}
    df = pd.DataFrame(data)
    for c in ["Position", "TimePoint", "Colony", "TrackNumber"]:
        if c in df:
            df[c] = df[c].fillna(0).astype(int)
    return df


def read_position_csv(path, columns=None):
    """Read one *_tracked.csv file from track_position.py (same role as the .mat files)."""
    df = pd.read_csv(path, usecols=lambda c: columns is None or c in columns)
    for c in ["Position", "TimePoint", "Colony", "TrackNumber"]:
        if c in df:
            df[c] = df[c].fillna(0).astype(int)
    return df
