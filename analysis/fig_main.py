"""Main figures (pages 3-4) and percentage summaries (pages 13-14).

Final_Figure_Mouse / Final_Figure_Human (one page each, legend included):
  a  the question;  b1 population level with published thresholds (LysoBrite, TMRM);  b2 other markers with the
  cut-offs of the original analysis (<=1.1x / >=1.8x);  c main result (ratio): % of cells dividing asymmetrically by
  how they were born;  d the same comparison for every measure;  e lineages asymmetric in all divisions vs random
Part1_Fig3_mouse_summary_percent / Part2_Fig3_human_summary_percent: percentage pairs per measure and marker.

Usage (inside the results folder, after steps 1-2):  python fig_main.py <figure folder>
"""
import sys, os, textwrap, warnings; warnings.filterwarnings('ignore'); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle, FancyArrowPatch
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
OUT = sys.argv[1] if len(sys.argv) > 1 else 'figures'; os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({'font.family': 'Liberation Sans', 'font.size': 7, 'axes.titlesize': 7.5, 'axes.labelsize': 7, 'xtick.labelsize': 6.3,
                     'ytick.labelsize': 6.3, 'legend.fontsize': 6.3, 'axes.linewidth': 0.6, 'axes.spines.top': False, 'axes.spines.right': False, 'pdf.fonttype': 42})
SYM, INT, ASY, INK, GREY, MORE, LESS = '#2B8CBE', '#B4B4B4', '#E34A33', '#1F2A36', '#9AA5B1', '#B5422E', '#F2A391'
MK = {'mouse': ['LysoBrite', 'Sca1', 'CD41', 'FcgR', 'CD105', 'CD201', 'CD71'], 'human': ['LysoBrite', 'TMRM', 'CD71', 'CD201', 'CD33', 'Inta3', 'CD34']}
LAB = {'FcgR': 'FcγR', 'Inta3': 'CD49c'}; lab = lambda m: LAB.get(m, m)
MEAS = ['log2ratio', 'pctdiff', 'absdiff', 'absdiff_adj']; MSHORT = ['Ratio', '% diff.', 'Abs.\ndiff.', 'Abs. diff.\nadjusted']
def letter(ax, l, x=-0.12, y=1.03): ax.text(x, y, l, transform=ax.transAxes, fontsize=10, fontweight='bold', va='bottom')
def save(fig, name): fig.savefig(f'{OUT}/{name}.pdf'); fig.savefig(f'{OUT}/{name}.png', dpi=300); plt.close(fig)

R = pd.read_csv('birth_vs_division.csv'); L3 = pd.read_csv('lineage_enrichment_3div_mouse.csv'); L2 = pd.read_csv('lineage_enrichment_2div.csv')
cells, div, pairs = pd.read_pickle('dataset.pkl')

def born_table(sp):
    """% dividing asymmetrically: born symmetric vs born asymmetric (pooled), per measure and marker."""
    rows = []
    for (me, mk), g in R[R.Species == sp].groupby(['Measure', 'Marker']):
        g = g.set_index('group')
        if not {'sym', 'asym_more', 'asym_less', 'log2OR'} <= set(g.index): continue
        n_a = g.loc['asym_more', 'n'] + g.loc['asym_less', 'n']
        p_a = (g.loc['asym_more', 'value'] * g.loc['asym_more', 'n'] + g.loc['asym_less', 'value'] * g.loc['asym_less', 'n']) / n_a
        o = g.loc['log2OR']
        rows.append(dict(Measure=me, Marker=mk, born_sym=100 * g.loc['sym', 'value'], born_asym=100 * p_a,
                         diff=100 * (p_a - g.loc['sym', 'value']), sig=(o.lo > 0) or (o.hi < 0)))
    return pd.DataFrame(rows)

def lineage_table(sp):
    L = L3 if sp == 'mouse' else L2[L2.Species == 'human']
    rows = []
    for (me, mk), g in L.groupby(['Measure', 'Marker']):
        g = g.set_index('stat')
        A, N = g.loc['all_asym'], g.loc['never_asym']
        rows.append(dict(Measure=me, Marker=mk, n=int(A.n_lineages), all_obs=100 * A.observed, all_exp=100 * A.expected, all_sig=A.p < 0.05,
                         sw_obs=100 * (1 - A.observed - N.observed), sw_exp=100 * (1 - A.expected - N.expected), sw_sig=g.loc['switches', 'p'] < 0.05))
    return pd.DataFrame(rows)

CMAP = LinearSegmentedColormap.from_list('rb', [SYM, '#FFFFFF', ASY])
def heat(ax, M, T, S, mk, title, vmax, show_y=True):
    ax.imshow(M, cmap=CMAP, norm=TwoSlopeNorm(0, -vmax, vmax), aspect='auto')
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            if T[i][j]: ax.text(j, i, T[i][j] + ('*' if S[i, j] else ''), ha='center', va='center', fontsize=5.4,
                                color='white' if abs(M[i, j]) > 0.65 * vmax else INK)
    ax.set_xticks(range(4)); ax.set_xticklabels(MSHORT, fontsize=5.8); ax.set_yticks(range(len(mk)))
    ax.set_yticklabels([lab(m) for m in mk] if show_y else []); ax.tick_params(length=0)
    [s.set_visible(False) for s in ax.spines.values()]; ax.set_title(title, fontsize=6.6)

def summary_mats(sp):
    mk = [m for m in MK[sp] if m in R[R.Species == sp].Marker.unique()]
    B, Lt = born_table(sp), lineage_table(sp)
    out = []
    for kind in ['born', 'all', 'sw']:
        M = np.full((len(mk), 4), np.nan); T = [[''] * 4 for _ in mk]; S = np.zeros((len(mk), 4), bool)
        for i, m in enumerate(mk):
            for j, me in enumerate(MEAS):
                if kind == 'born':
                    r = B[(B.Marker == m) & (B.Measure == me)]
                    if len(r): r = r.iloc[0]; M[i, j] = r['diff']; T[i][j] = f'{r.born_asym:.0f} vs {r.born_sym:.0f}'; S[i, j] = r.sig
                else:
                    r = Lt[(Lt.Marker == m) & (Lt.Measure == me)]
                    if len(r):
                        r = r.iloc[0]; o, e, sg = (r.all_obs, r.all_exp, r.all_sig) if kind == 'all' else (r.sw_obs, r.sw_exp, r.sw_sig)
                        M[i, j] = o - e; T[i][j] = f'{o:.0f} vs {e:.0f}'; S[i, j] = sg
        out.append((M, T, S))
    return mk, out

# ======================================================= summary figures with percentages
def summary(sp, name):
    mk, ((Mb, Tb, Sb), (Ma, Ta, Sa), (Ms, Ts, Ss)) = summary_mats(sp)
    nd = 3 if sp == 'mouse' else 2
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 3.1)); fig.subplots_adjust(left=0.11, right=0.99, top=0.78, bottom=0.2, wspace=0.12)
    heat(axes[0], Mb, Tb, Sb, mk, '% of cells dividing asymmetrically:\nborn asymmetric vs born symmetric', 20)
    heat(axes[1], Ma, Ta, Sa, mk, ('% of lineages asymmetric in all\n3 divisions' if nd == 3 else '% of lineages asymmetric in\nboth divisions') + ': observed vs random', 6, False)
    heat(axes[2], Ms, Ts, Ss, mk, f'% of lineages that switch between\nsym. and asym.: observed vs random', 12, False)
    for k, ax in enumerate(axes): fig.text(ax.get_position().x0 - (0.1 if k == 0 else 0.02), ax.get_position().y1 + 0.1, 'abc'[k], fontsize=10, fontweight='bold')
    fig.suptitle(f'{"Mouse" if sp == "mouse" else "Human"} HSC: summary per measure and marker', fontsize=8, fontweight='bold', y=0.985)
    fig.text(0.5, 0.04, 'Numbers: percentages (first vs second). Colour: difference (red = first higher, blue = first lower, white = no difference). '
             '* significant (95% CI or permutation p < 0.05). Chance = 33% in a.', ha='center', fontsize=5.8, color='#4A5560')
    save(fig, name)

# ======================================================= final one-page figures
PUB = {('human', 'LysoBrite'): (1.5, 3.0), ('human', 'TMRM'): (1.1, 1.5), ('mouse', 'LysoBrite'): (1.2, 1.5)}
def pub_classes(d):
    out = []
    for sp, mk, l in zip(d.Species, d.Marker, d.log2ratio):
        s_, a_ = PUB[(sp, mk)]; r = 2 ** l
        out.append('sym' if r < s_ else 'asym' if r > a_ else 'int')
    return np.array(out)

def mini_cartoon(ax):
    ax.axis('off'); ax.set_xlim(0, 12); ax.set_ylim(0, 5.6); ax.set_aspect('equal')
    def cell(x, y, a, r=0.38):
        g = 1 - a; ax.add_patch(Circle((x, y), r, facecolor=(0.15 + 0.85 * g, 0.62 + 0.38 * g, 0.36 + 0.64 * g), edgecolor=INK, lw=0.6, zorder=4))
    def div(x, y, a0, a1, a2, w, dy, kind, r=0.38):
        ax.plot([x, x, x - w, x - w], [y, y - dy * .45, y - dy * .45, y - dy], color=GREY, lw=0.8); ax.plot([x, x + w, x + w], [y - dy * .45, y - dy * .45, y - dy], color=GREY, lw=0.8)
        ax.add_patch(Circle((x, y - dy * .45), 0.16, color=ASY if kind == 'A' else SYM, zorder=3)); cell(x, y, a0, r); cell(x - w, y - dy, a1, r); cell(x + w, y - dy, a2, r)
    ax.text(2.3, 5.25, 'mother divides', ha='center', fontsize=6.3)
    div(2.3, 4.6, .6, .9, .15, 1.2, 2.2, 'A')
    ax.text(1.1, 1.85, 'received\nmore', ha='center', va='top', fontsize=6, color=MORE); ax.text(3.5, 1.85, 'received\nless', ha='center', va='top', fontsize=6, color='#C06A55')
    ax.add_patch(FancyArrowPatch((4.4, 2.6), (5.5, 2.6), arrowstyle='-|>', mutation_scale=9, color=INK, lw=0.9))
    ax.text(8.6, 5.25, 'when a daughter divides …', ha='center', fontsize=6.3)
    div(7.1, 4.3, .9, 1.0, .4, 0.75, 1.9, 'A', r=0.32); ax.text(7.1, 1.55, 'asymmetric again\n= passed on', ha='center', va='top', fontsize=5.8, color=ASY)
    ax.text(8.6, 3.0, 'or', fontsize=6.3, ha='center', va='center')
    div(10.1, 4.3, .9, .85, .85, 0.75, 1.9, 'S', r=0.32); ax.text(10.1, 1.55, 'symmetric or asym.\nat random', ha='center', va='top', fontsize=5.8, color=SYM)
    ax.set_title('Do daughters pass on what they received?', fontsize=7.2, fontweight='bold', loc='left')

def _pop_bars(ax, d, mk, sp, gens):
    w = 0.8 / len(gens)
    for i, m in enumerate(mk):
        for j, g in enumerate(gens):
            h = d[(d.Marker == m) & (d.DaughterGen == g)]
            if len(h) < 20: continue
            fr = h.cls.value_counts(normalize=True); x = i + (j - (len(gens) - 1) / 2) * w; b = 0
            for k, col in zip(['sym', 'int', 'asym'], [SYM, INT, ASY]):
                ax.bar(x, 100 * fr.get(k, 0), bottom=b, width=w * 0.9, color=col, edgecolor='white', linewidth=0.3); b += 100 * fr.get(k, 0)
            ax.text(x, 101.5, str(g), ha='center', fontsize=4.8, color='#4A5560')
    ax.set_xticks(range(len(mk))); ax.set_xticklabels([lab(m) for m in mk], rotation=30, ha='right'); ax.set_ylim(0, 108)
    ax.set_xlim(-0.6, len(mk) - 0.4)

def population(axes, sp):
    """b1: markers with published thresholds (Loeffler et al. 2019, 2022); b2: other markers with the cut-offs of the
    original analysis (<=1.1x symmetric, >=1.8x asymmetric), for description only."""
    ax1, ax2 = axes
    d = div[div.Species == sp].dropna(subset=['log2ratio']).copy(); r = 2 ** d.log2ratio
    gens = sorted(g for g in d.DaughterGen.unique() if (d.DaughterGen == g).sum() >= 300)
    pmk = [m for m in ['LysoBrite', 'TMRM'] if (sp, m) in PUB]
    dp = d[d.Marker.isin(pmk)].copy()
    lo = np.array([PUB[(sp, m)][0] for m in dp.Marker]); hi = np.array([PUB[(sp, m)][1] for m in dp.Marker]); rp = 2 ** dp.log2ratio.values
    dp['cls'] = np.where(rp < lo, 'sym', np.where(rp > hi, 'asym', 'int'))
    _pop_bars(ax1, dp, pmk, sp, gens); ax1.set_ylabel('% of divisions')
    ax1.set_xticklabels([f'{m}\n<{PUB[(sp, m)][0]:g}×, >{PUB[(sp, m)][1]:g}×' for m in pmk], rotation=0, ha='center', fontsize=5.2)
    ax1.set_title('b1  Published\n      thresholds', fontsize=6.6, loc='left')
    omk = [m for m in MK[sp] if m in d.Marker.unique() and m not in pmk]
    do = d[d.Marker.isin(omk)].copy(); ro = 2 ** do.log2ratio.values
    do['cls'] = np.where(ro <= 1.1, 'sym', np.where(ro >= 1.8, 'asym', 'int'))
    _pop_bars(ax2, do, omk, sp, gens); ax2.set_yticklabels([])
    ax2.set_title('b2  Other markers: no published threshold,\n      original cut-offs (≤1.1× / ≥1.8×)', fontsize=6.6, loc='left')
    ax2.legend(handles=[Rectangle((0, 0), 1, 1, color=c) for c in [SYM, INT, ASY]], labels=['symmetric', 'intermediate', 'asymmetric'],
               frameon=False, ncol=3, fontsize=5.5, loc='upper center', bbox_to_anchor=(0.5, -0.3), handlelength=1)

def born_bars(ax, sp):
    q = R[(R.Species == sp) & (R.Measure == 'log2ratio')]; mk = [m for m in MK[sp] if m in q.Marker.unique()]
    groups = [('sym', SYM, 'born symmetric'), ('asym_more', MORE, 'born asymmetric, received more'), ('asym_less', LESS, 'born asymmetric, received less')]
    for i, m in enumerate(mk):
        g = q[q.Marker == m].set_index('group')
        for k, (gr, col, _) in enumerate(groups):
            if gr not in g.index: continue
            v = g.loc[gr]; x = i + (k - 1) * 0.27
            ax.bar(x, 100 * v.value, width=0.25, color=col)
            ax.errorbar(x, 100 * v.value, yerr=[[100 * (v.value - v.lo)], [100 * (v.hi - v.value)]], color=INK, lw=0.5, capsize=0.8)
    ax.axhline(33.3, color=INK, ls='--', lw=0.6); ax.text(len(mk) - 0.45, 34.5, 'chance', fontsize=5.5, ha='right', color='#4A5560')
    ax.set_xticks(range(len(mk))); ax.set_xticklabels([lab(m) for m in mk]); ax.set_ylim(0, 75); ax.set_ylabel('% dividing asymmetrically')
    ax.set_title('Main result (ratio): how often cells divide asymmetrically, by how they were born', fontsize=7, loc='left')
    ax.legend(handles=[Rectangle((0, 0), 1, 1, color=c) for _, c, _ in groups], labels=[l for _, _, l in groups], frameon=False, ncol=3, fontsize=5.8, loc='upper left')

def lineage_bars(ax, sp):
    Lt = lineage_table(sp); Lt = Lt[Lt.Measure == 'log2ratio']; mk = [m for m in MK[sp] if m in Lt.Marker.unique()]
    nd = 3 if sp == 'mouse' else 2
    for i, m in enumerate(mk):
        r = Lt[Lt.Marker == m].iloc[0]
        ax.bar(i - 0.18, r.all_obs, width=0.34, color=ASY); ax.bar(i + 0.18, r.all_exp, width=0.34, color='white', edgecolor=INK, lw=0.6)
        if r.all_sig: ax.text(i, max(r.all_obs, r.all_exp) + 0.4, '*', ha='center', fontsize=7)
    ax.set_xticks(range(len(mk))); ax.set_xticklabels([lab(m) for m in mk], rotation=30, ha='right'); ax.set_ylabel('% of lineages')
    ax.set_ylim(0, max(Lt.all_obs.max(), Lt.all_exp.max()) * 1.35)
    ax.set_title('Lineages asymmetric in all 3 divisions (ratio)' if nd == 3 else 'Lineages asymmetric in both divisions (ratio)', fontsize=7, loc='left')
    ax.legend(handles=[Rectangle((0, 0), 1, 1, color=ASY), Rectangle((0, 0), 1, 1, facecolor='white', edgecolor=INK)], labels=['observed', 'expected if random'],
              frameon=False, fontsize=5.8, loc='upper right')

def final(sp, name, legend_text):
    fig = plt.figure(figsize=(7.2, 9.4))
    gs = fig.add_gridspec(4, 6, height_ratios=[1.05, 1.1, 1.3, 1.25], hspace=0.6, wspace=1.4, left=0.09, right=0.98, top=0.96, bottom=0.02)
    ax = fig.add_subplot(gs[0, 0:3]); mini_cartoon(ax); fig.text(0.005, ax.get_position().y1 + 0.012, 'a', fontsize=10, fontweight='bold')
    sub = gs[0, 3:6].subgridspec(1, 2, width_ratios=[1.1, 4] if sp == 'mouse' else [2.8, 4.4], wspace=0.12)
    ax1 = fig.add_subplot(sub[0]); ax2 = fig.add_subplot(sub[1]); population((ax1, ax2), sp)
    fig.text(ax1.get_position().x0 - 0.05, ax1.get_position().y1 + 0.012, 'b', fontsize=10, fontweight='bold')
    ax = fig.add_subplot(gs[1, :]); born_bars(ax, sp); fig.text(0.005, ax.get_position().y1 + 0.012, 'c', fontsize=10, fontweight='bold')
    mk, ((Mb, Tb, Sb), (Ma, Ta, Sa), _) = summary_mats(sp)
    ax = fig.add_subplot(gs[2, 0:3]); heat(ax, Mb, Tb, Sb, mk, '% dividing asymmetrically, born asym. vs sym.,\nfor each measure (colour: difference)', 20); fig.text(0.005, ax.get_position().y1 + 0.025, 'd', fontsize=10, fontweight='bold')
    ax = fig.add_subplot(gs[2, 3:6]); lineage_bars(ax, sp); letter(ax, 'e', -0.2, 1.1)
    ax = fig.add_subplot(gs[3, :]); ax.axis('off')
    ax.text(0, 1.0, legend_text[0], va='top', fontsize=7, fontweight='bold', transform=ax.transAxes)
    ax.text(0, 0.9, '\n'.join(textwrap.fill(p, 150) for p in legend_text[1:]), va='top', fontsize=6.2, linespacing=1.35, transform=ax.transAxes)
    save(fig, name)

def stats_sentence(sp):
    B = born_table(sp); r = B[B.Measure == 'log2ratio']
    return r

MOUSE_LEGEND = [
 "Mouse HSC: daughters partially pass on asymmetric inheritance",
 "(a) Question. A mother cell divides symmetrically or asymmetrically; one daughter receives more and the other less of a marker. When these daughters divide themselves, "
 "is their division asymmetric again (inheritance is passed on) or symmetric/asymmetric at random? "
 "(b) Population level: share of symmetric, intermediate and asymmetric divisions per daughter generation (1–3, numbers above bars). b1: LysoBrite with the published mouse thresholds (symmetric <1.2-fold, asymmetric >1.5-fold; Loeffler et al. 2019). b2: markers without a published threshold, shown for description with the cut-offs of the original analysis (≤1.1-fold / ≥1.8-fold); not directly comparable with b1. "
 "(c) Main result: percentage of cells that divide asymmetrically themselves, for cells born in a symmetric division (blue) or an asymmetric division in which they received more (dark red) "
 "or less (light red) of the marker. Asymmetry measured as the ratio between sisters; classes are thirds of all divisions of the same experiment, condition, marker and generation, so chance = 33% (dashed line). "
 "Error bars: 95% confidence intervals from resampling colonies. For all seven markers, cells born asymmetrically divide asymmetrically more often than cells born symmetrically. "
 "(d) The same comparison for every measure of asymmetry: ratio, % difference, absolute difference and absolute difference adjusted for the amount of marker. Numbers: % dividing asymmetrically when born "
 "asymmetric vs born symmetric; colour: difference (red = more often after asymmetric birth); * 95% CI of the odds ratio excludes 1. The absolute difference gives larger differences because it partly reflects "
 "how much marker a cell has, which is itself inherited. (e) Lineages of three consecutive divisions (e.g. divisions of cells 1, 2 and 4) that are asymmetric in all three: observed vs expected if divisions were random "
 "(classes shuffled within experiment × condition × marker × generation; * permutation p < 0.05).",
 "Data: curated lineage trees (generations 0–3) of 2 experiments (181117DL14, 181219DL13), 4 conditions each; marker values at the first measurement after birth."]
HUMAN_LEGEND = [
 "Human HSC: no consistent evidence that daughters pass on asymmetric inheritance",
 "(a) Question, as in the mouse figure. (b) Population level per daughter generation (1–2). b1: LysoBrite and TMRM with the published thresholds (LysoBrite <1.5-fold / >3-fold, TMRM <1.1-fold / >1.5-fold; Loeffler et al. 2022). b2: markers without a published threshold, shown for description with the cut-offs of the original analysis (≤1.1-fold / ≥1.8-fold); not directly comparable with b1. "
 "(c) Percentage of cells that divide asymmetrically themselves, for cells born symmetric (blue) or asymmetric, receiving more (dark red) or less (light red); ratio measure; classes are thirds of all divisions of "
 "the same experiment, condition, marker and generation, so chance = 33% (dashed). Error bars: 95% CI from resampling colonies. Differences are small and mostly within the confidence intervals. "
 "(d) The same comparison for all measures of asymmetry; numbers: % dividing asymmetrically when born asymmetric vs symmetric; colour: difference; * 95% CI of the odds ratio excludes 1. "
 "Effects appear with the absolute difference but disappear after adjusting for the amount of marker, indicating that they reflect inherited marker levels rather than inherited asymmetry. "
 "(e) Lineages asymmetric in both of two consecutive divisions (first → second): observed vs expected if random (* permutation p < 0.05).",
 "Data: curated lineage trees (generations 0–2; generation-2 cells rarely divide before the movie ends; HSC conditions) of 4 experiments (191125FS13, 191207FS14, 191216FS13, 191216FS14). Only one generation step is curated, so statistical power is "
 "lower than in mouse (≈100–1,300 mother–daughter pairs per marker)."]

summary('mouse', 'Part1_Fig3_mouse_summary_percent'); summary('human', 'Part2_Fig3_human_summary_percent')
final('mouse', 'Final_Figure_Mouse', MOUSE_LEGEND); final('human', 'Final_Figure_Human', HUMAN_LEGEND)
