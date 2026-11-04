"""Figures requested after the first presentation (pages 2 and 9-12).

  A_tracking_depth_and_n            example curated tree per species and n numbers
  B_pctdiff_distributions_by_birth  % difference without thresholds, by how the cell was born
  C_changes_over_generations        marker level, absolute and % difference per generation
  D_lineage_pattern_counts          how often each lineage pattern (A->A->A, ...) occurs vs random
  E_sisters_cousins_by_founder      do sisters and cousins divide alike, depending on the founder division?

Usage (inside the results folder, after steps 1 and 3):
  python fig_supplementary.py <figure folder> fig_depth fig_pctdiff fig_generations fig_patterns fig_sisters
"""
import sys, os, warnings; warnings.filterwarnings('ignore'); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
OUT = sys.argv[1] if len(sys.argv) > 1 else 'figures'; os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({'font.family': 'Liberation Sans', 'font.size': 7, 'axes.titlesize': 7.5, 'axes.labelsize': 7, 'xtick.labelsize': 6.3,
                     'ytick.labelsize': 6.3, 'legend.fontsize': 6.2, 'axes.linewidth': 0.6, 'axes.spines.top': False, 'axes.spines.right': False, 'pdf.fonttype': 42})
SYM, INT, ASY, INK, GREY, MORE, LESS = '#2B8CBE', '#B4B4B4', '#E34A33', '#1F2A36', '#9AA5B1', '#B5422E', '#F2A391'
MK = {'mouse': ['LysoBrite', 'Sca1', 'CD41', 'FcgR', 'CD105', 'CD201', 'CD71'], 'human': ['LysoBrite', 'TMRM', 'CD71', 'CD201', 'CD33', 'Inta3', 'CD34']}
LAB = {'FcgR': 'FcγR', 'Inta3': 'CD49c'}; lab = lambda m: LAB.get(m, m)
MCOL = ['#1F2A36', '#E34A33', '#2B8CBE', '#6A51A3', '#31A354', '#E6A400', '#8C564B']
def letter(fig, ax, l, dx=0.06, dy=0.015): p = ax.get_position(); fig.text(max(p.x0 - dx, 0.005), p.y1 + dy, l, fontsize=10, fontweight='bold')
def save(fig, name): fig.savefig(f'{OUT}/{name}.pdf'); fig.savefig(f'{OUT}/{name}.png', dpi=300); plt.close(fig)
cells, div, pairs = pd.read_pickle('dataset.pkl')
NN = pd.read_csv('n_numbers.csv').set_index('Species')

# ======================================================= A: tracking depth, example trees, n
def draw_tree(ax, sp, exp, pos, col, marker='LysoBrite'):
    t = cells[(cells.Experiment == exp) & (cells.Position == pos) & (cells.Colony == col)].drop_duplicates('TrackNumber').set_index('TrackNumber')
    d = div[(div.Experiment == exp) & (div.Position == pos) & (div.Colony == col) & (div.Marker == marker)].set_index('Mother')
    tracks = set(t.index); xp = {}; c = [0]
    def place(n):
        kids = [k for k in (2 * n, 2 * n + 1) if k in tracks]
        if not kids: xp[n] = c[0]; c[0] += 1
        else:
            for k in kids: place(k)
            xp[n] = np.mean([xp[k] for k in kids])
    place(1)
    for n in xp:
        r = t.loc[n]; ax.plot([xp[n]] * 2, [r.FirstTP, r.LastTP], color='#7B8794', lw=0.9)
        kids = [k for k in (2 * n, 2 * n + 1) if k in xp]
        if kids:
            ax.plot([min(xp[k] for k in kids), max(xp[k] for k in kids)], [r.LastTP] * 2, color='#7B8794', lw=0.9)
            if n in d.index:
                k = d.loc[n, 'log2ratio_class']; ax.scatter(xp[n], r.LastTP, s=16, color={'sym': SYM, 'int': INT, 'asym': ASY}.get(k, 'k'), zorder=3, edgecolor='white', lw=0.4)
        if r.Fate == 'died': ax.text(xp[n], r.LastTP + 4, '×', ha='center', va='top', fontsize=6, color=INK)
    gens = t.groupby('Generation').FirstTP.min(); xm = max(xp.values())
    for g, y in gens.items(): ax.text(xm + 0.9, y, f'gen {g}', fontsize=5.8, color='#4A5560', va='center', ha='left')
    ax.set_xlim(-0.8, xm + 0.9 + 0.12 * (xm + 2)); ax.set_ylim(t.LastTP.max() + 12, t.FirstTP.min() - 8); ax.set_xticks([])
    ax.spines['bottom'].set_visible(False); ax.set_ylabel('time (frames)')
    return t

def fig_depth():
    fig = plt.figure(figsize=(7.2, 6.4)); gs = fig.add_gridspec(2, 2, height_ratios=[1.5, 1], hspace=0.3, wspace=0.35, left=0.1, right=0.98, top=0.93, bottom=0.07)
    # example trees: the most complete curated tree per species
    for k, sp in enumerate(['mouse', 'human']):
        c = cells[cells.Species == sp].copy()
        c['tree'] = list(zip(c.Experiment, c.Position, c.Colony))
        g = c.groupby('tree').agg(n=('TrackNumber', 'nunique'), maxg=('Generation', 'max'), div=('Fate', lambda s: (s == 'divided').sum()))
        if sp == 'human':
            g2 = g[(g.maxg == 2) & (g.n >= 7)]
            best = g2.sort_values('div', ascending=False).index[0]
        else:
            best = g[g.maxg == 3].sort_values('div', ascending=False).index[0]
        ax = fig.add_subplot(gs[0, k]); draw_tree(ax, sp, *best); letter(fig, ax, 'ab'[k])
        if k == 0: ax0 = ax
        ax.set_title(f'{"Mouse" if sp == "mouse" else "Human"}: example curated tree ({best[0]}, position {best[1]}, colony {best[2]})', fontsize=6.6)
    ax0.legend(handles=[plt.Line2D([], [], marker='o', ls='', color=c, ms=4) for c in [SYM, INT, ASY]],
                                             labels=['symmetric', 'intermediate', 'asymmetric (LysoBrite)'], frameon=False, loc='lower left', fontsize=5.8)
    # n table
    ax = fig.add_subplot(gs[1, :]); ax.axis('off'); letter(fig, ax, 'c', 0.08, -0.02)
    rows = [('Experiments', 'experiments'), ('Conditions analysed', 'conditions'), ('Trees (founder cells) tracked', 'trees_tracked'), ('Trees with ≥ 1 division', 'trees_with_division')]
    rows += [(f'Cells in generation {g}  (of which divided)', f'gen{g}') for g in range(4)]
    tab = []
    for name, key in rows:
        line = [name]
        for sp in ['mouse', 'human']:
            r = NN.loc[sp]
            if key.startswith('gen'):
                g = int(key[3:]); line.append(f"{int(r[f'cells_gen{g}']):,}  ({int(r[f'divided_gen{g}']):,})" if r[f'cells_gen{g}'] > 0 else '–')
            else: line.append(f'{int(r[key]):,}')
        tab.append(line)
    tb = ax.table(cellText=tab, colLabels=['', 'Mouse HSC', 'Human HSC'], colWidths=[0.42, 0.22, 0.22], loc='upper center', cellLoc='center')
    tb.auto_set_font_size(False); tb.set_fontsize(6.2); tb.scale(1, 1.0)
    for (i, j), cl in tb.get_celld().items():
        cl.set_edgecolor('#D9DDE2')
        if i == 0: cl.set_facecolor('#1F2A36'); cl.get_text().set_color('white'); cl.get_text().set_fontweight('bold')
        if j == 0: cl.get_text().set_ha('left')
    fig.text(0.06, 0.012, 'Human: generation-2 cells are tracked (3,345) but almost never divide within the movie (28): most reach the end of the movie (1,175) or their tracking ends.\n'
            'So human trees are limited by movie length, not by curation. Mouse: generations 0–3 used (three rounds of division), as curated.',
            fontsize=5.8, color='#4A5560', va='bottom')
    fig.suptitle('How deep were the experiments tracked?', fontsize=8.5, fontweight='bold', y=0.985)
    save(fig, 'A_tracking_depth_and_n')

# ======================================================= B: % difference distributions by birth
def fig_pctdiff():
    raw = pd.read_pickle('pctdiff_raw.pkl'); st = pd.read_csv('pctdiff_by_birth.csv')
    fig, axes = plt.subplots(2, 1, figsize=(7.2, 5.4)); fig.subplots_adjust(left=0.08, right=0.99, top=0.89, bottom=0.08, hspace=0.45)
    order = [('sym', SYM, 'born symmetric'), ('int', INT, 'born intermediate'), ('asym_more', MORE, 'born asym., received more'), ('asym_less', LESS, 'born asym., received less')]
    for ax, sp in zip(axes, ['mouse', 'human']):
        mk = [m for m in MK[sp] if m in raw[raw.c_Species == sp].Marker.unique()]
        for i, m in enumerate(mk):
            for k, (b, col, _) in enumerate(order):
                v = raw[(raw.c_Species == sp) & (raw.Marker == m) & (raw.born == b)].c_pctdiff.clip(upper=200)
                if len(v) < 10: continue
                x = i + (k - 1.5) * 0.2
                vp = ax.violinplot([v], [x], widths=0.19, showextrema=False)
                for bd in vp['bodies']: bd.set_facecolor(col); bd.set_alpha(0.75); bd.set_edgecolor('none')
                ax.plot([x - 0.07, x + 0.07], [v.median()] * 2, color=INK, lw=1)
            r = st[(st.Species == sp) & (st.Marker == m)].iloc[0]
            sig = r.lo > 0 or r.hi < 0
            ax.text(i, 207, f'{r.shift_pct_points:+.0f}' + ('*' if sig else ''), ha='center', fontsize=6.3, color=ASY if (sig and r.shift_pct_points > 0) else '#7A8591', fontweight='bold' if sig else None)
        ax.set_xticks(range(len(mk))); ax.set_xticklabels([lab(m) for m in mk]); ax.set_ylim(0, 218); ax.set_ylabel("% difference in the cell's\nown division")
        ax.set_title(f'{"Mouse" if sp == "mouse" else "Human"}', loc='left', fontsize=7.5, fontweight='bold')
    axes[0].legend(handles=[Rectangle((0, 0), 1, 1, color=c) for _, c, _ in order], labels=[l for _, _, l in order], ncol=4, frameon=False, loc='upper center', bbox_to_anchor=(0.5, 1.17))
    fig.suptitle('Distribution of % difference, by how the cell was born (no threshold)', fontsize=8.5, fontweight='bold', y=0.995)
    fig.text(0.5, 0.01, 'Violins: distribution; line: median. Number above each marker: how much more unequal (in percentile points) the divisions of asymmetrically born cells are than those of\n'
             'symmetrically born cells (0 = no shift; * 95% CI excludes 0). Birth class = thirds of the mother division (% difference) within experiment × condition × marker × generation.', ha='center', fontsize=5.8, color='#4A5560')
    save(fig, 'B_pctdiff_distributions_by_birth')

# ======================================================= C: changes over generations
def fig_generations():
    g = pd.read_csv('generation_trends.csv')
    fig, axes = plt.subplots(2, 3, figsize=(7.2, 4.8)); fig.subplots_adjust(left=0.09, right=0.84, top=0.9, bottom=0.1, wspace=0.42, hspace=0.5)
    for r, sp in enumerate(['mouse', 'human']):
        mk = MK[sp]
        for c, (col, title) in enumerate([('total', 'Marker per division\n(both daughters, median)'), ('absdiff', 'Absolute difference\nbetween sisters (median)'), ('pctdiff', '% difference\nbetween sisters (median)')]):
            ax = axes[r, c]
            for i, m in enumerate(mk):
                h = g[(g.Species == sp) & (g.Marker == m) & (g.n >= 50)].sort_values('DaughterGen')
                if h.empty: continue
                y = h[col] / h[col].iloc[0] if col != 'pctdiff' else h[col]
                ax.plot(h.DaughterGen, y, marker='o', ms=2.5, lw=1, color=MCOL[i], label=lab(m))
            ax.set_xticks(sorted(g[(g.Species == sp) & (g.n >= 50)].DaughterGen.unique())); ax.set_xlabel('daughter generation')
            if col != 'pctdiff': ax.axhline(1, color=GREY, lw=0.5, ls='--'); ax.set_ylabel('relative to generation 1')
            else: ax.set_ylabel('% difference')
            if r == 0: ax.set_title(title, fontsize=7)
        axes[r, 2].legend(frameon=False, fontsize=5.8, loc='center left', bbox_to_anchor=(1.02, 0.5), title='Mouse' if sp == 'mouse' else 'Human', title_fontsize=6.5)
        fig.text(0.01, 0.72 - r * 0.43, 'Mouse' if sp == 'mouse' else 'Human', rotation=90, fontsize=8, fontweight='bold', va='center')
    fig.suptitle('Across generations, the absolute difference follows the marker level; the % difference much less', fontsize=8, fontweight='bold', y=0.99)
    save(fig, 'C_changes_over_generations')

# ======================================================= D: lineage pattern counts
def fig_patterns():
    pm = pd.read_csv('patterns_mouse.csv'); ph = pd.read_csv('patterns_human.csv')
    fig = plt.figure(figsize=(7.2, 7.4)); gs = fig.add_gridspec(3, 2, height_ratios=[1, 1.15, 1], width_ratios=[1.6, 1], hspace=0.6, wspace=0.3, left=0.1, right=0.98, top=0.93, bottom=0.07)
    def bars(ax, d, title):
        pats = list(d.pattern); x = np.arange(len(pats))
        ax.bar(x - 0.2, d['count'], width=0.4, color=[ASY if set(p) == {'A'} else SYM if set(p) == {'S'} else '#6E7B88' for p in pats], label='observed')
        ax.bar(x + 0.2, d.expected_count, width=0.4, color='white', edgecolor=INK, lw=0.6, label='expected if random')
        for xi, r in zip(x, d.itertuples()):
            ax.text(xi, max(r.count, r.expected_count) * 1.02, f'{r.count}' + ('*' if r.p < 0.05 else ''), ha='center', fontsize=5.3, rotation=90 if len(pats) > 8 else 0, va='bottom')
        ax.set_xticks(x); ax.set_xticklabels(['→'.join(p) for p in pats], rotation=90 if len(pats) > 8 else 0, fontsize=6); ax.set_ylabel('number of lineages'); ax.set_title(title, loc='left', fontsize=7)
        ax.set_ylim(0, max(d['count'].max(), d.expected_count.max()) * 1.4)
    order = lambda d: d.assign(nA=d.pattern.str.count('A')).sort_values(['nA', 'pattern'], ascending=[False, True])
    d3 = order(pm[(pm.Measure == 'log2ratio') & (pm.Marker == 'LysoBrite') & (pm.length == 3)])
    ax = fig.add_subplot(gs[0, 0]); bars(ax, d3, f'Mouse, LysoBrite: lineages of 3 divisions (n = {int(d3.n_lineages.iloc[0]):,})'); ax.legend(handles=[Rectangle((0, 0), 1, 1, color=ASY), Rectangle((0, 0), 1, 1, color='#6E7B88'), Rectangle((0, 0), 1, 1, color=SYM),
                       Rectangle((0, 0), 1, 1, facecolor='white', edgecolor=INK, lw=0.6)],
              labels=['observed: asymmetric in every division', 'observed: mixed (switches)', 'observed: symmetric in every division', 'expected if random'],
              frameon=False, loc='upper left', ncol=2, fontsize=5.8); letter(fig, ax, 'a')
    dh = order(ph[(ph.Measure == 'log2ratio') & (ph.Marker == 'LysoBrite')])
    ax = fig.add_subplot(gs[0, 1]); bars(ax, dh, f'Human, LysoBrite: 2 divisions (n = {int(dh.n_lineages.iloc[0]):,})'); letter(fig, ax, 'b', 0.04)
    # heatmaps: all markers, observed % vs expected %
    def heat(ax, d, mk, title):
        pats = list(order(d[d.Marker == mk[0]]).pattern)
        M = np.full((len(mk), len(pats)), np.nan); T = [[''] * len(pats) for _ in mk]; S = np.zeros_like(M, bool)
        for i, m in enumerate(mk):
            for j, p in enumerate(pats):
                r = d[(d.Marker == m) & (d.pattern == p)]
                if len(r): r = r.iloc[0]; M[i, j] = r.pct - r.expected_pct; T[i][j] = f'{r.pct:.0f}\n{r.expected_pct:.0f}'; S[i, j] = r.p < 0.05
        ax.imshow(M, cmap=LinearSegmentedColormap.from_list('rb', [SYM, '#FFFFFF', ASY]), norm=TwoSlopeNorm(0, -6, 6), aspect='auto')
        for i in range(len(mk)):
            for j in range(len(pats)):
                if T[i][j]: ax.text(j, i, T[i][j] + ('*' if S[i, j] else ''), ha='center', va='center', fontsize=4.9, linespacing=0.95)
        ax.set_xticks(range(len(pats))); ax.set_xticklabels(['→'.join(p) for p in pats], fontsize=6); ax.set_yticks(range(len(mk))); ax.set_yticklabels([lab(m) for m in mk])
        ax.tick_params(length=0); [s.set_visible(False) for s in ax.spines.values()]; ax.set_title(title, loc='left', fontsize=7)
    ax = fig.add_subplot(gs[1, :]); heat(ax, pm[(pm.Measure == 'log2ratio') & (pm.length == 3)], MK['mouse'],
                                         'Mouse, all markers (ratio): % of lineages with each pattern, observed (top) vs expected if random (bottom)'); letter(fig, ax, 'c')
    ax = fig.add_subplot(gs[2, 0]); heat(ax, ph[ph.Measure == 'log2ratio'], MK['human'], 'Human, all markers (ratio): observed (top) vs expected (bottom), %'); letter(fig, ax, 'd')
    ax = fig.add_subplot(gs[2, 1]); ax.axis('off')
    ax.text(0, 0.95, 'How to read', fontsize=7, fontweight='bold', va='top')
    ax.text(0, 0.8, 'Each lineage (divisions of cells 1, 2, 4, …) is written\nas one letter per division, from the founder down:\nA = more unequal than the median division,\nS = more equal than the median.\n\nA→A→A = asymmetric in all three divisions,\nS→S→S = symmetric in all three.\n\nRed: seen more often than expected if\ndivisions were random; blue: less often.\n* p < 0.05 (permutation within experiment ×\ncondition × marker × generation).',
            fontsize=6, va='top', color='#4A5560', linespacing=1.3)
    fig.suptitle('How often does each lineage behaviour occur?', fontsize=8.5, fontweight='bold', y=0.985)
    save(fig, 'D_lineage_pattern_counts')

# ======================================================= E: sisters and cousins by founder division
def fig_sisters():
    s = pd.read_csv('sisters_cousins.csv'); s = s[(s.Measure == 'log2ratio') & (s.n_pairs >= 30)]
    rels = [('sisters (cells 2 vs 3)', 'Sisters\n(daughters of the founder)'), ('sisters (gen 2: 4v5, 6v7)', 'Sisters, next generation\n(granddaughters of the founder)'), ('cousins (4/5 vs 6/7)', 'Cousins\n(granddaughters, different mothers)')]
    fig, axes = plt.subplots(2, 3, figsize=(7.2, 5.0), sharey=True); fig.subplots_adjust(left=0.1, right=0.98, top=0.86, bottom=0.14, wspace=0.12, hspace=0.45)
    for r, sp in enumerate(['mouse', 'human']):
        mk = MK[sp]
        for c, (rel, title) in enumerate(rels):
            ax = axes[r, c]; q = s[(s.Species == sp) & (s.relation == rel)]
            for i, m in enumerate(mk):
                for k, (f, col) in enumerate([('sym', SYM), ('int', INT), ('asym', ASY)]):
                    v = q[(q.Marker == m) & (q.founder == f)]
                    if not len(v): continue
                    v = v.iloc[0]; x = i + (k - 1) * 0.22
                    ax.errorbar(x, v.pct_same, yerr=[[v.pct_same - v.lo], [v.hi - v.pct_same]], fmt='o', ms=2.8, color=col, lw=0.6, capsize=0.8)
            ax.axhline(33.3, color=INK, ls='--', lw=0.6); ax.set_xticks(range(len(mk))); ax.set_xticklabels([lab(m) for m in mk], rotation=45, ha='right', fontsize=5.8)
            ax.set_ylim(0, 75)
            if r == 0: ax.set_title(title, fontsize=6.8)
            if c == 0: ax.set_ylabel(f'{"Mouse" if sp == "mouse" else "Human"}\n% of pairs dividing\nthe same way')
            if q.empty: ax.text(0.5, 0.5, 'not enough curated\ngenerations', transform=ax.transAxes, ha='center', fontsize=6.5, color='#7A8591')
    fig.legend(handles=[plt.Line2D([], [], marker='o', ls='', color=c, ms=4) for c in [SYM, INT, ASY]] + [plt.Line2D([], [], color=INK, ls='--', lw=0.6)],
               labels=['founder division symmetric', 'intermediate', 'asymmetric', 'chance (33%)'], ncol=4, frameon=False, loc='lower center', fontsize=6.2)
    fig.suptitle('Do sisters and cousins divide alike, depending on how the founder divided?', fontsize=8.5, fontweight='bold', y=0.985)
    fig.text(0.5, 0.925, 'Same way = both divisions in the same class (symmetric / intermediate / asymmetric thirds, ratio). Error bars: 95% CI (colony resampling).', ha='center', fontsize=6, color='#4A5560')
    save(fig, 'E_sisters_cousins_by_founder')

if __name__ == '__main__':
    for f in sys.argv[2:]: globals()[f]()
