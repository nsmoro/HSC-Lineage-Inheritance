"""Robustness check (page 17): mouse results with generation 4 included.

Compares the main result (% of cells dividing asymmetrically, born asymmetric vs born symmetric; ratio, thirds)
between the main analysis (mouse generations 0-3) and a second run with MAX_GEN_MOUSE=4.

Usage (inside the results folder):  python fig_generation4_check.py <figure folder> <results folder of the generation-4 run>
Output: G5_check_generation4.pdf / .png
"""
import os, sys, textwrap, warnings
warnings.filterwarnings("ignore")
import pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "Liberation Sans"})
OUT, GEN4 = sys.argv[1], sys.argv[2]
LAB = {"FcgR": "FcγR"}
ORDER = ["LysoBrite", "Sca1", "CD41", "FcgR", "CD105", "CD201", "CD71"]


def born_asym_vs_sym(path):
    """Per mouse marker: (% dividing asymmetrically if born asymmetric, if born symmetric, significant?)."""
    r = pd.read_csv(path)
    r = r[(r.Species == "mouse") & (r.Measure == "log2ratio")]
    out = {}
    for mk, g in r.groupby("Marker"):
        g = g.set_index("group")
        n = g.loc["asym_more", "n"] + g.loc["asym_less", "n"]
        pa = (g.loc["asym_more", "value"] * g.loc["asym_more", "n"] + g.loc["asym_less", "value"] * g.loc["asym_less", "n"]) / n
        sig = g.loc["log2OR", "lo"] > 0 or g.loc["log2OR", "hi"] < 0
        out[mk] = f"{100 * pa:.0f}% vs {100 * g.loc['sym', 'value']:.0f}%" + ("*" if sig else "")
    return out


a = born_asym_vs_sym("birth_vs_division.csv")
b = born_asym_vs_sym(os.path.join(GEN4, "birth_vs_division.csv"))
rows = [[LAB.get(m, m), a[m], b[m]] for m in ORDER if m in a and m in b]
fig = plt.figure(figsize=(7.2, 4.2))
fig.text(0.05, 0.93, "Check: mouse results with generation 4 included", fontsize=11, fontweight="bold")
ax = fig.add_axes([0.08, 0.3, 0.84, 0.55]); ax.axis("off")
tb = ax.table(cellText=rows, colLabels=["Marker", "Generations 0–3 (main analysis)", "Generations 0–4"], loc="center", cellLoc="center")
tb.auto_set_font_size(False); tb.set_fontsize(8); tb.scale(1, 1.3)
for (i, j), c in tb.get_celld().items():
    c.set_edgecolor("#D9DDE2")
    if i == 0:
        c.set_facecolor("#1F2A36"); c.get_text().set_color("white"); c.get_text().set_fontweight("bold")
fig.text(0.05, 0.2, textwrap.fill("% of cells dividing asymmetrically, born asymmetric vs born symmetric (ratio, thirds; chance 33%); "
         "* 95% CI of the odds ratio excludes 1. The main analysis uses generations 0–3, where the curation signatures (explicit cell "
         "fates, both daughters tracked) are strongest.", 125), fontsize=6.8, va="top", linespacing=1.4)
fig.savefig(f"{OUT}/G5_check_generation4.pdf"); fig.savefig(f"{OUT}/G5_check_generation4.png", dpi=300)
print("ok")
