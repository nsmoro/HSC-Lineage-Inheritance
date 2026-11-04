"""Figure: the question (page 1). Follows the draft by D. Loeffler.

Panels: (1) lineage tree with TrackNumbers, (2) distribution of sister asymmetry over all divisions (example: mouse
LysoBrite) and the idea of following the daughters, (3) every lineage as a row of A/S, (4) analysis split by species,
measure and marker, (5) how divisions are classified (thirds; published thresholds for population plots).

Usage (inside the results folder, after step 1):  python fig_cartoon.py <figure folder>
Output: Fig0_cartoon_DL_draft.pdf / .png
"""
import sys, os, warnings; warnings.filterwarnings('ignore')
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyBboxPatch, FancyArrowPatch, Rectangle
plt.rcParams.update({'font.family': 'Liberation Sans', 'font.size': 7, 'pdf.fonttype': 42})
SYM, INT, ASY, INK, GREY = '#2B8CBE', '#B4B4B4', '#E34A33', '#1F2A36', '#9AA5B1'
OUT = sys.argv[1] if len(sys.argv) > 1 else 'figures'; os.makedirs(OUT, exist_ok=True)

fig = plt.figure(figsize=(7.2, 7.6))
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 14.4); ax.set_ylim(-3.2, 12); ax.axis('off'); ax.set_aspect('equal')
def arrow(p0, p1, rad=0.0, lw=0.9, col=INK, ms=8):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle='-|>', mutation_scale=ms, color=col, lw=lw, connectionstyle=f'arc3,rad={rad}'))
def panel(x, y, w, h, num, title):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.04,rounding_size=0.2', facecolor='#F6F7F9', edgecolor='#D9DDE2', lw=0.8))
    ax.text(x + 0.2, y + h - 0.22, num, fontsize=11, fontweight='bold', color=ASY, va='top')
    ax.text(x + 0.62, y + h - 0.25, title, fontsize=7.8, fontweight='bold', va='top', color=INK)

# ---------------- 1: lineage tree over time
panel(0.15, 5.0, 5.6, 6.85, '1', 'Track every cell → lineage tree')
cls = {1: 'A', 2: 'A', 3: 'S', 4: 'A', 5: 'A', 6: 'S', 7: 'A'}
xs = {}; top, dy = 10.2, 1.45
def place(n, x, w):
    xs[n] = x
    if n < 8: place(2 * n, x - w, w / 2); place(2 * n + 1, x + w, w / 2)
place(1, 3.3, 1.12)
for n, x in xs.items():
    g = int(np.floor(np.log2(n))); y = top - g * dy
    if n < 8:
        ym = y - dy * 0.5
        ax.plot([x, x], [y - 0.2, ym], color=GREY, lw=0.8); ax.plot([xs[2 * n], xs[2 * n + 1]], [ym, ym], color=GREY, lw=0.8)
        for k in (2 * n, 2 * n + 1): ax.plot([xs[k], xs[k]], [ym, y - dy + 0.2], color=GREY, lw=0.8)
        ax.add_patch(Circle((x, ym), 0.11, color=ASY if cls[n] == 'A' else SYM, zorder=4))
    ax.add_patch(Rectangle((x - 0.22, y - 0.19), 0.44, 0.38, facecolor='white', edgecolor=INK, lw=0.7, zorder=5))
    ax.text(x, y, str(n), ha='center', va='center', fontsize=5.6 if n > 9 else 6.3, zorder=6)
for g in range(4): ax.text(0.32, top - g * dy, f'gen {g}', fontsize=6.5, va='center', color='#4A5560')
arrow((5.5, 10.4), (5.5, 6.4), lw=0.8); ax.text(5.62, 8.4, 'time', rotation=-90, fontsize=6.5, va='center', color='#4A5560')
ax.text(2.95, 5.2, 'boxes: TrackNumber (daughters of n = 2n, 2n+1)\ndots: division, red asymmetric / blue symmetric', ha='center', fontsize=5.8, color='#4A5560')

# ---------------- 2: distribution -> pick sym / asym divisions
panel(5.95, 5.0, 8.3, 6.85, '2', 'Score every division, then follow the daughters')
cells, div, pairs = pd.read_pickle('dataset.pkl')
v = div[(div.Species == 'mouse') & (div.Marker == 'LysoBrite')].log2ratio.dropna(); v = v[v <= 3]
lo, hi = np.percentile(v, [100 / 3, 200 / 3])
h = fig.add_axes([0.475, (8.46 + 3.2) / 15.2, 0.27, 2.28 / 15.2])
n, b, patches = h.hist(v, bins=60, color=INT, edgecolor='none')
for p, left in zip(patches, b[:-1]):
    p.set_facecolor(SYM if left < lo else ASY if left >= hi else INT)
h.set_xlabel('asymmetry between sisters\n(e.g. ratio, log2)', fontsize=6); h.set_ylabel('# divisions', fontsize=6)
h.tick_params(labelsize=5.5); h.spines['top'].set_visible(False); h.spines['right'].set_visible(False)
ymax = n.max()
h.annotate('', xy=(0, ymax * 1.08), xytext=(lo, ymax * 1.08), arrowprops=dict(arrowstyle='-', color=SYM, lw=1.6), annotation_clip=False)
h.text(lo / 2, ymax * 1.14, 'Sym', color=SYM, fontsize=7, fontweight='bold', ha='center')
h.annotate('', xy=(hi, ymax * 0.55), xytext=(3, ymax * 0.55), arrowprops=dict(arrowstyle='-', color=ASY, lw=1.6))
h.text((hi + 3) / 2, ymax * 0.6, 'Asym', color=ASY, fontsize=7, fontweight='bold', ha='center'); h.set_xlim(0, 3)
h.set_title('all divisions (example: mouse, LysoBrite)', fontsize=6, color='#4A5560')
# side text
ax.text(11.4, 10.6, 'most even third\n→ symmetric (S)', fontsize=6.2, color=SYM, va='top')
ax.text(11.4, 9.7, 'middle third\n→ intermediate', fontsize=6.2, color='#7A7A7A', va='top')
ax.text(11.4, 8.8, 'most unequal third\n→ asymmetric (A)', fontsize=6.2, color=ASY, va='top')
# small trees: follow the daughters
def small_tree(x0, y0, kind, label, col):
    ax.add_patch(Circle((x0, y0), 0.16, facecolor='#7FBF96', edgecolor=INK, lw=0.6, zorder=4))
    ax.plot([x0, x0], [y0 - 0.16, y0 - 0.45], color=GREY, lw=0.8); ax.plot([x0 - 0.55, x0 + 0.55], [y0 - 0.45, y0 - 0.45], color=GREY, lw=0.8)
    ax.add_patch(Circle((x0, y0 - 0.45), 0.09, color=col, zorder=5))
    a1, a2 = (.9, .15) if kind == 'A' else (.55, .55)
    for dx, a in [(-0.55, a1), (0.55, a2)]:
        g_ = 1 - a
        ax.plot([x0 + dx, x0 + dx], [y0 - 0.45, y0 - 0.75], color=GREY, lw=0.8)
        ax.add_patch(Circle((x0 + dx, y0 - 0.9), 0.15, facecolor=(0.15 + 0.85 * g_, 0.62 + 0.38 * g_, 0.36 + 0.64 * g_), edgecolor=INK, lw=0.6, zorder=4))
        arrow((x0 + dx, y0 - 1.1), (x0 + dx * 1.5, y0 - 1.55), rad=0.25 if dx < 0 else -0.25, lw=0.8)
        ax.text(x0 + dx * 1.6, y0 - 1.85, '?', ha='center', fontsize=9, fontweight='bold')
    ax.text(x0, y0 + 0.35, label, ha='center', fontsize=6.3, color=col, fontweight='bold')
small_tree(7.3, 7.15, 'S', 'a symmetric division', SYM); small_tree(11.6, 7.15, 'A', 'an asymmetric division', ASY)
ax.text(9.45, 5.18, 'then follow the daughters: how do they divide next?', ha='center', fontsize=6.2, color='#4A5560')

# ---------------- 3: lineage table
panel(0.15, 0.15, 9.2, 4.65, '3', 'Write every lineage as a row: one division type per generation')
rows = [('1-2-4-8', 'AAA'), ('1-2-4-9', 'AAA'), ('1-2-5-10', 'AAS'), ('1-2-5-11', 'AAS'), ('1-3-6-12', 'ASA'), ('1-3-6-13', 'ASA'), ('1-3-7-14', 'ASA'), ('1-3-7-15', 'ASA')]
pat = {'1-2-4-8': 'AAA', '1-2-4-9': 'AAA', '1-2-5-10': 'AAA', '1-2-5-11': 'AAA', '1-3-6-12': 'ASS', '1-3-6-13': 'ASS', '1-3-7-14': 'ASA', '1-3-7-15': 'ASA'}
# use the classes of the tree in panel 1
def lin_pat(path):
    c = [int(x) for x in path.split('-')][:-1]
    return ''.join(cls[k] for k in c)
x0, y0, cw, rh = 2.6, 3.65, 0.62, 0.36
for j, g in enumerate(['gen 1', 'gen 2', 'gen 3']): ax.text(x0 + j * cw + cw / 2, y0 + 0.28, g, ha='center', fontsize=6.3, fontweight='bold')
ax.text(x0 + 3 * cw + 0.1, y0 + 0.28, '(gen 4 …)', fontsize=5.8, color='#4A5560')
paths = ['1-2-4-8', '1-2-4-9', '1-2-5-10', '1-2-5-11', '1-3-6-12', '1-3-6-13', '1-3-7-14', '1-3-7-15']
for i, pth in enumerate(paths):
    y = y0 - (i + 0.5) * rh
    ax.text(x0 - 0.12, y, f'lineage {pth}', ha='right', va='center', fontsize=5.8)
    for j, c in enumerate(lin_pat(pth)):
        ax.add_patch(Rectangle((x0 + j * cw, y - rh / 2), cw, rh, facecolor=ASY if c == 'A' else SYM, edgecolor='white', lw=1))
        ax.text(x0 + j * cw + cw / 2, y, c, ha='center', va='center', color='white', fontsize=6.5, fontweight='bold')
ax.plot([0.5, 0.5], [y0 - 8 * rh + 0.05, y0 - 0.05], color=INK, lw=1.2); ax.text(0.4, y0 - 4 * rh, '1 tree', rotation=90, ha='right', va='center', fontsize=6.3, fontweight='bold')
ax.text(x0 + 3 * cw + 0.35, y0 - 1.0, 'Question in one line:\ndoes a row that starts with A\nstay A (A A A), or does it switch\n(A S A, A A S) as often as\nexpected by chance?', fontsize=6.3, va='top', color=INK)
ax.text(x0 + 3 * cw + 0.35, y0 - 2.75, 'gen = generation of the dividing cell\n(divisions of cells 1, 2, 4 for lineage 1-2-4-8)', fontsize=5.6, va='top', color='#4A5560')
# ---------------- 4: divide by
panel(9.55, 0.15, 4.7, 4.65, '4', 'Done separately for')
items = [('1', 'Mouse', '2 experiments, 7 markers'), ('2', 'Human', '4 experiments, HSC, 7 markers'),
         ('3', 'Measurements', 'ratio · % difference ·\nabsolute difference (± level)'), ('4', 'Markers', 'each marker on its own')]
for i, (k, t, d) in enumerate(items):
    y = 3.75 - i * 0.92
    ax.add_patch(Circle((10.05, y), 0.2, color=INK)); ax.text(10.05, y, k, color='white', ha='center', va='center', fontsize=6.5, fontweight='bold')
    ax.text(10.4, y + 0.1, t, fontsize=7, fontweight='bold', va='center'); ax.text(10.4, y - 0.22, d, fontsize=5.8, va='top', color='#4A5560')

# ---------------- 5: how symmetric / asymmetric are defined (thirds)
panel(0.15, -3.1, 14.1, 3.05, '5', 'How a division is called symmetric or asymmetric: split into thirds')
folds = [1.02, 1.05, 1.08, 1.2, 1.3, 1.5, 1.9, 2.4, 3.1]
x0, yc, step = 1.0, -1.55, 0.95
for i, f in enumerate(folds):
    x = x0 + i * step
    t = np.tanh(np.log(f))                                # how unequal the two sisters are
    for dx, amt in [(-0.17, 0.55 + 0.4 * t), (0.17, 0.55 - 0.45 * t)]:
        g_ = 1 - amt
        ax.add_patch(Circle((x + dx, yc), 0.15, facecolor=(0.15 + 0.85 * g_, 0.62 + 0.38 * g_, 0.36 + 0.64 * g_), edgecolor=INK, lw=0.5, zorder=4))
    ax.text(x, yc - 0.38, f'{f:g}×', ha='center', fontsize=6, color=INK)
ax.text(x0 - 0.45, yc, 'sister\npairs', ha='right', va='center', fontsize=5.8, color='#4A5560')
arrow((x0 - 0.2, yc - 0.75), (x0 + 8 * step + 0.2, yc - 0.75), lw=0.7)
ax.text(x0 + 4 * step, yc - 0.98, 'sorted from most equal to most unequal sisters (number = fold difference between the two sisters)', ha='center', fontsize=5.8, color='#4A5560')
for k, (lab_, col) in enumerate([('most even third\n= symmetric', SYM), ('middle third\n= intermediate', '#7A7A7A'), ('most unequal third\n= asymmetric', ASY)]):
    xa, xb = x0 + 3 * k * step - 0.38, x0 + (3 * k + 2) * step + 0.38
    ax.plot([xa, xa, xb, xb], [yc + 0.3, yc + 0.38, yc + 0.38, yc + 0.3], color=col, lw=1.2)
    ax.text((xa + xb) / 2, yc + 0.45, lab_, ha='center', va='bottom', fontsize=6.1, color=col, fontweight='bold')
ax.text(10.0, -0.7, 'Three groups of equal size', fontsize=6.6, fontweight='bold', va='top')
ax.text(10.0, -1.0, '• done separately for each experiment,\n  condition, marker and generation\n• if divisions were random, 1 in 3 would be\n  asymmetric in any group → chance = 33%\n• works the same for every measure\n  (ratio, % difference, absolute difference)',
        fontsize=5.8, va='top', color='#4A5560', linespacing=1.25)
ax.text(x0 - 0.7, -2.95, 'Population plots use the published thresholds (mouse LysoBrite <1.2× / >1.5×; human LysoBrite <1.5× / >3×, TMRM <1.1× / >1.5×); lineage patterns (A/S) use halves.',
        fontsize=5.5, color='#4A5560')
fig.savefig(f'{OUT}/Fig0_cartoon_DL_draft.pdf'); fig.savefig(f'{OUT}/Fig0_cartoon_DL_draft.png', dpi=300); print('ok')
