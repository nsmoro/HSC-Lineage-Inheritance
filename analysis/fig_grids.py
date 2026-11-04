"""Detail figures (pages 15-16): every measure (rows) x every marker (columns).

Bars: % of cells dividing asymmetrically by how they were born (symmetric, intermediate, asymmetric received more,
asymmetric received less), 95% CI by colony resampling, chance = 33%. Label: odds ratio born asymmetric vs symmetric.

Usage (inside the results folder, after step 2):  python fig_grids.py <figure folder>
"""
import sys, os, warnings; warnings.filterwarnings('ignore'); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyBboxPatch, Rectangle, FancyArrowPatch
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm, LogNorm
OUT = sys.argv[1] if len(sys.argv) > 1 else 'figures'; os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({'font.family': 'Liberation Sans', 'font.size': 7, 'axes.titlesize': 7.5, 'axes.labelsize': 7, 'xtick.labelsize': 6.3,
                     'ytick.labelsize': 6.3, 'legend.fontsize': 6.3, 'axes.linewidth': 0.6, 'axes.spines.top': False, 'axes.spines.right': False, 'pdf.fonttype': 42})
SYM, INT, ASY, INK, GREY, MORE, LESS = '#2B8CBE', '#B4B4B4', '#E34A33', '#1F2A36', '#9AA5B1', '#B5422E', '#F2A391'
MK = {'mouse': ['LysoBrite', 'Sca1', 'CD41', 'FcgR', 'CD105', 'CD201', 'CD71'], 'human': ['LysoBrite', 'TMRM', 'CD71', 'CD201', 'CD33', 'Inta3', 'CD34']}
LAB = {'FcgR': 'FcγR', 'Inta3': 'CD49c'}; lab = lambda m: LAB.get(m, m)
MEAS = ['log2ratio', 'pctdiff', 'absdiff', 'absdiff_adj']
MLAB = {'log2ratio': 'Ratio', 'pctdiff': '% difference', 'absdiff': 'Absolute difference', 'absdiff_adj': 'Abs. difference,\nadjusted for level'}
def letter(ax, l, x=-0.12, y=1.03): ax.text(x, y, l, transform=ax.transAxes, fontsize=10, fontweight='bold', va='bottom')
def save(fig, name): fig.savefig(f'{OUT}/{name}.pdf'); fig.savefig(f'{OUT}/{name}.png', dpi=300); plt.close(fig)
def cell(ax, x, y, amount, r=0.16, z=4):
    g = 1 - amount; ax.add_patch(Circle((x, y), r, facecolor=(0.15 + 0.85 * g, 0.62 + 0.38 * g, 0.36 + 0.64 * g), edgecolor=INK, lw=0.6, zorder=z))
def division(ax, x, y, amt, a1, a2, w, dy, kind, r=0.16, lw=0.8):
    ax.plot([x, x, x - w, x - w], [y, y - dy * 0.45, y - dy * 0.45, y - dy], color=GREY, lw=lw, zorder=1)
    ax.plot([x, x + w, x + w], [y - dy * 0.45, y - dy * 0.45, y - dy], color=GREY, lw=lw, zorder=1)
    ax.add_patch(Circle((x, y - dy * 0.45), r * 0.4, color=ASY if kind == 'A' else SYM, zorder=3))
    cell(ax, x, y, amt, r); cell(ax, x - w, y - dy, a1, r); cell(ax, x + w, y - dy, a2, r)
def arrow(ax, p0, p1, **kw): ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle='-|>', mutation_scale=8, color=INK, lw=0.9, **kw))

# =============================================================== per species: grid measure x marker
R = pd.read_csv('birth_vs_division.csv')
def grid(sp, name):
    mk = [m for m in MK[sp] if m in R[R.Species == sp].Marker.unique()]
    fig, axes = plt.subplots(4, len(mk), figsize=(7.2, 6.4), sharey=True, sharex=True)
    fig.subplots_adjust(left=0.11, right=0.99, top=0.9, bottom=0.11, wspace=0.12, hspace=0.25)
    groups = [('sym', SYM, 'born symmetric'), ('int', INT, 'born intermediate'), ('asym_more', MORE, 'born asymmetric, received more'), ('asym_less', LESS, 'born asymmetric, received less')]
    for i, me in enumerate(MEAS):
        for j, m in enumerate(mk):
            ax = axes[i, j]; q = R[(R.Species == sp) & (R.Measure == me) & (R.Marker == m)].set_index('group')
            for k, (g, col, _) in enumerate(groups):
                if g in q.index:
                    v = q.loc[g]; ax.bar(k, 100 * v.value, color=col, width=0.75)
                    ax.errorbar(k, 100 * v.value, yerr=[[100 * (v.value - v.lo)], [100 * (v.hi - v.value)]], color=INK, lw=0.5, capsize=0.8)
            ax.axhline(33.3, color=INK, ls='--', lw=0.5)
            if 'log2OR' in q.index:
                o = q.loc['log2OR']; sig = o.lo > 0 or o.hi < 0
                ax.text(0.5, 0.97, ('↑ ' if o.value > 0 else '↓ ') + f'{2 ** o.value:.1f}×' + ('*' if sig else ''), transform=ax.transAxes, ha='center', va='top', fontsize=5.6,
                        color=ASY if (sig and o.value > 0) else SYM if (sig and o.value < 0) else '#7A8591')
            ax.set_xticks([]); ax.set_ylim(0, 75)
            if i == 0: ax.set_title(lab(m), fontsize=7, fontweight='bold')
            if j == 0: ax.set_ylabel(MLAB[me], fontsize=6.5)
    fig.text(0.015, 0.5, '% of cells that divide asymmetrically themselves', rotation=90, va='center', fontsize=7)
    fig.legend(handles=[Rectangle((0, 0), 1, 1, color=c) for _, c, _ in groups] + [plt.Line2D([], [], color=INK, ls='--', lw=0.6)],
               labels=[l for _, _, l in groups] + ['chance (33%)'], loc='lower center', ncol=5, frameon=False, fontsize=6.2)
    fig.suptitle(f'{"Mouse" if sp == "mouse" else "Human"} HSC: do daughters pass on what they received?  (rows: measure, columns: marker)', fontsize=8, fontweight='bold', y=0.975)
    fig.text(0.5, 0.935, 'numbers: odds of dividing asymmetrically when born asymmetric vs born symmetric (↑ more, ↓ less often; * 95% CI excludes 1)', ha='center', fontsize=6, color='#4A5560')
    save(fig, name)
grid('mouse', 'Part1_Fig2_mouse_by_measure_and_marker')
grid('human', 'Part2_Fig2_human_by_measure_and_marker')

