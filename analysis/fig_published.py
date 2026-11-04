"""Figures for the analyses aligned with the published study (pages 5-8).

  F1  published thresholds: share of asymmetric / symmetric divisions (per experiment, by generation)
  F2  published thresholds: is lysosome / TMRM asymmetry passed on?
  F3  co-inheritance of each marker with lysosomes (LysoHigh daughter)
  F4  how LysoHigh and LysoLow daughters divide themselves
Each figure is saved as PNG (300 dpi) and as a vector page with title and legend (<name>_page.pdf).

Usage (inside the results folder, after step 4):  python fig_published.py <figure folder>
"""
import sys
import warnings; warnings.filterwarnings('ignore')
import numpy as np, pandas as pd, matplotlib, os
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.family': 'Liberation Sans', 'font.size': 11, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.linewidth': 0.8})
SYM, INT, ASY, INK, GREY = '#2B6CB0', '#A7AFB8', '#D9534F', '#1E2832', '#8A96A3'
HI, LO = '#8C2D04', '#F5A65B'
OUT = sys.argv[1] if len(sys.argv) > 1 else 'figures'; os.makedirs(OUT, exist_ok=True)
LAB = {'Inta3': 'CD49c', 'FcgR': 'FcγR'}; lab = lambda m: LAB.get(m, m)
rng = np.random.default_rng(1)
import textwrap
PAGES = {
 'F1_published_thresholds_population': ('Published thresholds: how many divisions are asymmetric',
   ['Left: share of first divisions classified as asymmetric (red) or symmetric (blue) with the published thresholds (human LysoBrite >3× / <1.5× and TMRM >1.5× / <1.1×, Loeffler et al. 2022; mouse LysoBrite >1.5× / <1.2×, Loeffler et al. 2019). Dots: independent experiments; bars: mean ± SEM. Diamond: value reported in the paper for human LysoBrite (27.1 ± 3.6%). Right: the same classes by daughter generation (mean of experiments).',
    'Human LysoBrite: 13.1 ± 1.9% asymmetric (4 experiments, 1,538 divisions), lower than published; likely reasons are the measurement window after division (paper: up to 60 min) and the set of experiments analysed.']),
 'F2_published_thresholds_inheritance': ('Published thresholds: is lysosome asymmetry passed on?',
   ["Percentage of cells that divide asymmetrically themselves, by how they were born, with the published thresholds for both divisions. Dots: experiments; bars: mean ± SEM. Tests: two-sided Mann–Whitney on the fold change of the cell's own division (all cells); paired t-test across experiments (human, 4 experiments).",
    'Mouse LysoBrite: 44% vs 38% after asymmetric vs symmetric birth (p < 0.001). Human LysoBrite: 12% vs 11% (n.s.). Human TMRM: 28% vs 20% (p = 0.07 across experiments).']),
 'F3_lysosome_coinheritance': ('Markers are co-inherited with lysosomes; active mitochondria segregate to the other daughter',
   ['For every asymmetric lysosome division (published thresholds): percentage of divisions in which each marker is higher in the LysoHigh daughter (50% = random). Dots: experiments; bars: mean ± SEM.',
    'Human: CD71 98%, CD201 99%, CD49c 100%, CD33 87%, CD34 96%; TMRM 32% (68% go to the LysoLow daughter). Mouse markers: 77–95%. This reproduces the co-inheritance and anti-segregation reported by Loeffler et al. 2022.']),
 'F4_lysohigh_lysolow_daughters': ('LysoLow daughters divide more asymmetrically themselves',
   ['Cells born from a symmetric lysosome division (blue), or as the LysoHigh (dark brown) or LysoLow (orange) daughter of an asymmetric lysosome division: percentage whose own division is asymmetric for each marker (LysoBrite and TMRM: published thresholds; other markers: most unequal third).',
    'Human CD71: 41% vs 19% (p = 0.006; consistent in all 4 experiments, paired t p = 0.012). Mouse: LysoBrite 46% vs 42%, CD71 40% vs 33%, CD105 42% vs 35%, FcγR 39% vs 33% (all p < 0.05). For human LysoBrite and CD71 and mouse LysoBrite, CD71 and CD105 the effect remains after adjusting for marker level; for mouse FcγR it does not.'])}
def save(fig, n):
    fig.savefig(f'{OUT}/{n}.png', dpi=300, bbox_inches='tight', facecolor='white')
    if n in PAGES:
        t, leg = PAGES[n]
        fig.suptitle(t, x=0.01, ha='left', y=1.08, fontsize=13, fontweight='bold')
        fig.text(0.01, -0.12 if 'F4' not in n else -0.2, '\n'.join(textwrap.fill(p, 150) for p in leg), fontsize=9, va='top', linespacing=1.45)
        fig.savefig(f'{OUT}/{n}_page.pdf', bbox_inches='tight', pad_inches=0.4)
    plt.close(fig)
def dots(ax, x, vals, col, w=0.12, bw=0.22):
    vals = np.asarray(vals, float); vals = vals[~np.isnan(vals)]
    if not len(vals): return
    ax.scatter(x + rng.uniform(-w, w, len(vals)), vals, s=26, color=col, edgecolor='white', lw=0.6, zorder=3)
    m = vals.mean(); se = vals.std(ddof=1) / np.sqrt(len(vals)) if len(vals) > 1 else 0
    ax.plot([x - bw, x + bw], [m, m], color=INK, lw=2, zorder=4); ax.errorbar(x, m, yerr=se, color=INK, lw=1.2, capsize=4, zorder=4)

# ---- F1: how many divisions are asymmetric with published thresholds (per experiment)
pop = pd.read_csv('pub_population.csv')
fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), gridspec_kw={'width_ratios': [1, 1.4]})
ax = axes[0]; items = [('human', 'LysoBrite', '>3× / <1.5×'), ('human', 'TMRM', '>1.5× / <1.1×'), ('mouse', 'LysoBrite', '>1.5× / <1.2×')]
for i, (sp, mk, t) in enumerate(items):
    q = pop[(pop.Species == sp) & (pop.Marker == mk) & (pop.DaughterGen == 1)]
    dots(ax, i - 0.18, q.asym, ASY); dots(ax, i + 0.18, q.sym, SYM)
ax.scatter([0 - 0.18], [27.1], marker='D', s=60, facecolor='white', edgecolor=ASY, lw=1.8, zorder=5)
ax.annotate('paper: 27.1 ± 3.6%', (-0.18, 27.1), xytext=(0.15, 36), fontsize=9.5, color=ASY, arrowprops=dict(arrowstyle='-', color=ASY, lw=0.8))
ax.set_xticks(range(3)); ax.set_xticklabels([f'{sp.capitalize()}\n{mk}\n{t}' for sp, mk, t in items], fontsize=9.5)
ax.set_ylabel('% of first divisions'); ax.set_ylim(0, 65)
ax.legend(handles=[plt.Line2D([], [], marker='o', ls='', color=ASY, label='asymmetric'), plt.Line2D([], [], marker='o', ls='', color=SYM, label='symmetric'),
                   plt.Line2D([], [], marker='D', ls='', mfc='white', mec=ASY, label='published value')], frameon=False, fontsize=9, loc='upper right')
ax.set_title('First division (dots = experiments, bar = mean ± SEM)', fontsize=10.5, loc='left')
ax = axes[1]
for j, (sp, mk, t) in enumerate(items):
    q = pop[(pop.Species == sp) & (pop.Marker == mk)].groupby('DaughterGen')[['sym', 'int', 'asym']].mean()
    q = q[pop[(pop.Species == sp) & (pop.Marker == mk)].groupby('DaughterGen').n.sum() >= 50]
    for k, (g, r) in enumerate(q.iterrows()):
        x = j * 5 + k; b = 0
        for c, col in [('sym', SYM), ('int', INT), ('asym', ASY)]:
            ax.bar(x, r[c], bottom=b, color=col, width=0.8); b += r[c]
        ax.text(x, -6, str(g), ha='center', fontsize=9)
    ax.text(j * 5 + (len(q) - 1) / 2, 106, f'{sp} {mk}', ha='center', fontsize=9.5, fontweight='bold')
ax.set_xticks([]); ax.set_ylim(-9, 114); ax.set_ylabel('% of divisions'); ax.text(-0.9, -6, 'gen', fontsize=9, ha='right')
ax.set_title('By daughter generation (mean of experiments)', fontsize=10.5, loc='left')
ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=c) for c in [SYM, INT, ASY]], labels=['symmetric', 'intermediate', 'asymmetric'], frameon=False, fontsize=9, ncol=3, loc='upper center', bbox_to_anchor=(0.5, -0.04))
save(fig, 'F1_published_thresholds_population')

# ---- F2: inheritance with published thresholds (per experiment)
inh = pd.read_csv('pub_inheritance_per_experiment.csv'); summ = pd.read_csv('pub_inheritance_summary.csv')
fig, axes = plt.subplots(1, 3, figsize=(11, 3.9), sharey=False)
for ax, (sp, mk, t) in zip(axes, items):
    q = inh[(inh.c_Species == sp) & (inh.Marker == mk)]
    for i, (b, col) in enumerate([('sym', SYM), ('int', INT), ('asym', ASY)]):
        dots(ax, i, q[q.born == b].pct_asym, col)
    s = summ[(summ.Species == sp) & (summ.Marker == mk)].iloc[0]
    ptxt = f'MWU p = {s.mwu_p:.2g}' + (f'\npaired t (experiments) p = {s.paired_t_p_experiments:.2g}' if not np.isnan(s.paired_t_p_experiments) else '\n(2 experiments: no test across experiments)')
    ax.text(0.02, 0.97, ptxt, transform=ax.transAxes, fontsize=8.5, va='top', color='#4A5560')
    ax.set_xticks(range(3)); ax.set_xticklabels(['born\nsymmetric', 'born\nintermediate', 'born\nasymmetric'], fontsize=9.5)
    ax.set_title(f'{sp.capitalize()} {mk} ({t})', fontsize=10.5); ax.set_ylim(0, max(70, q.pct_asym.max() * 1.25))
axes[0].set_ylabel('% of cells dividing\nasymmetrically themselves')
save(fig, 'F2_published_thresholds_inheritance')

# ---- F3: co-inheritance with lysosomes
co = pd.read_csv('lyso_coinheritance_per_experiment.csv')
fig, axes = plt.subplots(1, 2, figsize=(11, 3.9), gridspec_kw={'width_ratios': [6, 6]})
for ax, sp in zip(axes, ['human', 'mouse']):
    q = co[(co.Species == sp) & (co.L_class == 'asym')]
    mk = [m for m in (['CD71', 'CD201', 'Inta3', 'CD33', 'CD34', 'TMRM'] if sp == 'human' else ['Sca1', 'CD41', 'FcgR', 'CD105', 'CD201', 'CD71']) if m in q.Marker.unique()]
    for i, m in enumerate(mk):
        dots(ax, i, q[q.Marker == m].pct_co, '#5B8DB8' if m != 'TMRM' else HI)
    ax.axhline(50, color=GREY, ls='--', lw=1); ax.text(len(mk) - 0.5, 52, 'random', fontsize=9, color=GREY, ha='right')
    ax.set_xticks(range(len(mk))); ax.set_xticklabels([lab(m) for m in mk]); ax.set_ylim(0, 105)
    ax.set_title(f'{sp.capitalize()}: asymmetric lysosome divisions', fontsize=10.5, loc='left')
axes[0].set_ylabel('% of divisions where the marker\ngoes to the LysoHigh daughter')
save(fig, 'F3_lysosome_coinheritance')

# ---- F4: LysoHigh vs LysoLow daughters
lh = pd.read_csv('lyso_daughters_per_experiment.csv'); ls_ = pd.read_csv('lyso_daughters_summary.csv')
fig, axes = plt.subplots(1, 2, figsize=(11, 4.1), gridspec_kw={'width_ratios': [6, 7]})
fig.text(0.5, -0.13, 'Dots: experiments with ≥ 10 cells in the group; bar: mean ± SEM. * Mann–Whitney p < 0.05, LysoHigh vs LysoLow (all cells pooled).', ha='center', fontsize=9, color='#4A5560')
for ax, sp in zip(axes, ['human', 'mouse']):
    mk = ['LysoBrite', 'CD71', 'TMRM', 'CD201', 'Inta3', 'CD33'] if sp == 'human' else ['LysoBrite', 'CD71', 'CD105', 'CD201', 'CD41', 'FcgR', 'Sca1']
    q = lh[lh.c_Species == sp]
    for i, m in enumerate(mk):
        for k, (dname, col) in enumerate([('Lyso-symmetric birth', SYM), ('LysoHigh daughter', HI), ('LysoLow daughter', LO)]):
            v = q[(q.Marker == m) & (q.daughter == dname) & (q.n >= 10)].pct_asym
            dots(ax, i + (k - 1) * 0.27, v, col, w=0.04, bw=0.1)
        s = ls_[(ls_.Species == sp) & (ls_.Marker == m)]
        if len(s) and s.mwu_p_high_vs_low.iloc[0] < 0.05: ax.text(i, ax.get_ylim()[1] if False else 63, '*', ha='center', fontsize=14)
    ax.set_xticks(range(len(mk))); ax.set_xticklabels([lab(m) for m in mk]); ax.set_ylim(0, 68)
    ax.set_title(f'{sp.capitalize()}', fontsize=10.5, loc='left', fontweight='bold')
axes[0].set_ylabel('% of cells dividing\nasymmetrically themselves')
fig.legend(handles=[plt.Line2D([], [], marker='o', ls='', color=c) for c in [SYM, HI, LO]],
           labels=['born in a symmetric lysosome division', 'LysoHigh daughter of an asymmetric division', 'LysoLow daughter of an asymmetric division'],
           ncol=3, frameon=False, fontsize=9.5, loc='lower center', bbox_to_anchor=(0.5, -0.07))
save(fig, 'F4_lysohigh_lysolow_daughters')
print('ok')
