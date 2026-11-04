"""Step 4 - analyses aligned with the published study.

Conventions taken from Loeffler, Schneiter et al., Blood 2022 (human) and Loeffler et al., Nature 2019 (mouse):
  (1) marker-specific thresholds for symmetric / asymmetric divisions
        human LysoBrite <1.5x / >3x, human TMRM <1.1x / >1.5x, mouse LysoBrite <1.2x / >1.5x
  (4) daughters of an asymmetric lysosome division called LysoHigh / LysoLow; other markers read relative to them
  (5) one value per independent experiment (mean +- SEM); Mann-Whitney tests on pooled cells and paired t-tests
      across experiments (human, 4 experiments)

Usage (inside the results folder, after step 1):  python 04_published_conventions.py
Output:
  pub_population.csv, pub_first_division.csv                 division classes per experiment and generation
  pub_inheritance_per_experiment.csv, pub_inheritance_summary.csv   is lysosome / TMRM asymmetry passed on?
  lyso_coinheritance_per_experiment.csv, lyso_coinheritance_summary.csv   does each marker go to the LysoHigh daughter?
  lyso_daughters_per_experiment.csv, lyso_daughters_summary.csv     how LysoHigh / LysoLow daughters divide themselves
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import warnings; warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu, ttest_rel
cells, div, pairs = pd.read_pickle('dataset.pkl')
KEY = ['Experiment', 'Condition', 'Position', 'Colony']
# published thresholds: (symmetric below, asymmetric above), fold change between sisters
PUB = {('human', 'LysoBrite'): (1.5, 3.0), ('human', 'TMRM'): (1.1, 1.5), ('mouse', 'LysoBrite'): (1.2, 1.5)}

def pub_class(sp, mk, l2r):
    s, a = PUB[(sp, mk)]; r = 2 ** l2r
    return np.where(np.isnan(l2r), None, np.where(r < s, 'sym', np.where(r > a, 'asym', 'int')))

# ---------------------------------------------------------------- (1) published thresholds
d = div[[(s, m) in PUB for s, m in zip(div.Species, div.Marker)]].dropna(subset=['log2ratio']).copy()
d['pub'] = [pub_class(s, m, np.array([v]))[0] for s, m, v in zip(d.Species, d.Marker, d.log2ratio)]
pop = (d.groupby(['Species', 'Marker', 'Experiment', 'DaughterGen']).pub.value_counts(normalize=True).unstack(fill_value=0) * 100).reset_index()
pop_n = d.groupby(['Species', 'Marker', 'Experiment', 'DaughterGen']).size().rename('n').reset_index()
pop = pop.merge(pop_n); pop.to_csv('pub_population.csv', index=False)
# first-division summary across experiments (mean +- SEM), as in the paper (Fig 1D: 27.1 +- 3.6% > 3x)
f1 = pop[pop.DaughterGen == 1].groupby(['Species', 'Marker']).agg(asym_mean=('asym', 'mean'), asym_sem=('asym', 'sem'), sym_mean=('sym', 'mean'),
                                                                   sym_sem=('sym', 'sem'), n_exp=('asym', 'size'), n_div=('n', 'sum')).reset_index()
f1.to_csv('pub_first_division.csv', index=False)

# inheritance with published thresholds: how a cell divides, by how it was born
p = pairs[[(s, m) in PUB for s, m in zip(pairs.c_Species, pairs.Marker)]].dropna(subset=['p_log2ratio', 'c_log2ratio']).copy()
p['born'] = [pub_class(s, m, np.array([v]))[0] for s, m, v in zip(p.c_Species, p.Marker, p.p_log2ratio)]
p['own'] = [pub_class(s, m, np.array([v]))[0] for s, m, v in zip(p.c_Species, p.Marker, p.c_log2ratio)]
p['own_asym'] = (p.own == 'asym').astype(float)
inh = (p.groupby(['c_Species', 'Marker', 'Experiment', 'born']).agg(pct_asym=('own_asym', lambda x: 100 * x.mean()), n=('own_asym', 'size'),
                                                                  ratio_median=('c_log2ratio', lambda x: 2 ** x.median())).reset_index())
inh.to_csv('pub_inheritance_per_experiment.csv', index=False)
rows = []
for (sp, mk), g in p.groupby(['c_Species', 'Marker']):
    a = g[g.born == 'asym'].c_log2ratio; s = g[g.born == 'sym'].c_log2ratio
    mw = mannwhitneyu(a, s, alternative='two-sided') if len(a) > 5 and len(s) > 5 else None
    e = inh[(inh.c_Species == sp) & (inh.Marker == mk)].pivot(index='Experiment', columns='born', values='pct_asym')
    tt = ttest_rel(e['asym'], e['sym'], nan_policy='omit') if {'asym', 'sym'} <= set(e.columns) and len(e.dropna(subset=['asym','sym'])) >= 3 else None
    rows.append(dict(Species=sp, Marker=mk, n_born_asym=len(a), n_born_sym=len(s),
                     pct_asym_after_asym=100 * (g[g.born == 'asym'].own == 'asym').mean(), pct_asym_after_sym=100 * (g[g.born == 'sym'].own == 'asym').mean(),
                     median_fold_after_asym=2 ** a.median(), median_fold_after_sym=2 ** s.median(),
                     mwu_p=mw.pvalue if mw else np.nan, n_experiments=len(e), paired_t_p_experiments=tt.pvalue if tt else np.nan))
pd.DataFrame(rows).to_csv('pub_inheritance_summary.csv', index=False)

# ---------------------------------------------------------------- (4) LysoHigh / LysoLow reference
L = div[div.Marker == 'LysoBrite'][KEY + ['Mother', 'Species', 'log2ratio', 'high_daughter']].rename(columns={'log2ratio': 'L_l2r', 'high_daughter': 'L_high'}).dropna()
L['L_class'] = [pub_class(s, 'LysoBrite', np.array([v]))[0] for s, v in zip(L.Species, L.L_l2r)]
# 4a co-/anti-inheritance of each marker with lysosomes (division level)
m = div[div.Marker != 'LysoBrite'].dropna(subset=['log2ratio']).merge(L, on=KEY + ['Mother', 'Species'])
m['co'] = (m.high_daughter == m.L_high).astype(float)
co = (m.groupby(['Species', 'Marker', 'Experiment', 'L_class']).agg(pct_co=('co', lambda x: 100 * x.mean()), n=('co', 'size')).reset_index())
co.to_csv('lyso_coinheritance_per_experiment.csv', index=False)
co_sum = co.groupby(['Species', 'Marker', 'L_class']).agg(mean=('pct_co', 'mean'), sem=('pct_co', 'sem'), n_exp=('pct_co', 'size'), n=('n', 'sum')).reset_index()
co_sum.to_csv('lyso_coinheritance_summary.csv', index=False)
# 4b how LysoHigh vs LysoLow daughters divide themselves (each marker)
q = pairs.merge(L.rename(columns={'Mother': 'p_Mother', 'Species': 'c_Species'}), on=KEY + ['p_Mother', 'c_Species'])
q['daughter'] = ''
q.loc[q.L_class == 'sym', 'daughter'] = 'Lyso-symmetric birth'
q.loc[q.L_class == 'int', 'daughter'] = 'Lyso-intermediate birth'
hi = (q.c_Mother % 2 == 0) == (q.L_high == 2)
q.loc[(q.L_class == 'asym') & hi, 'daughter'] = 'LysoHigh daughter'
q.loc[(q.L_class == 'asym') & ~hi, 'daughter'] = 'LysoLow daughter'
# own division: LysoBrite/TMRM with published thresholds, other markers: top third (tertiles) as before
own = []
for r in q[['c_Species', 'Marker', 'c_log2ratio', 'c_log2ratio_class']].itertuples(index=False):
    if (r.c_Species, r.Marker) in PUB: own.append(pub_class(r.c_Species, r.Marker, np.array([r.c_log2ratio]))[0])
    else: own.append(r.c_log2ratio_class)
q['own'] = own; q = q.dropna(subset=['own'])
q['own_asym'] = (q.own == 'asym').astype(float)
lh = (q.groupby(['c_Species', 'Marker', 'Experiment', 'daughter']).agg(pct_asym=('own_asym', lambda x: 100 * x.mean()), n=('own_asym', 'size'),
                                                                     own_pct_rank=('c_log2ratio_pct', 'mean')).reset_index())
lh.to_csv('lyso_daughters_per_experiment.csv', index=False)
rows = []
for (sp, mk), g in q.groupby(['c_Species', 'Marker']):
    a = g[g.daughter == 'LysoHigh daughter']; b = g[g.daughter == 'LysoLow daughter']; s = g[g.daughter == 'Lyso-symmetric birth']
    mw = mannwhitneyu(a.c_log2ratio_pct.dropna(), b.c_log2ratio_pct.dropna(), alternative='two-sided') if len(a) > 5 and len(b) > 5 else None
    rows.append(dict(Species=sp, Marker=mk, n_high=len(a), n_low=len(b), n_sym=len(s), pct_asym_high=100 * a.own_asym.mean(), pct_asym_low=100 * b.own_asym.mean(),
                     pct_asym_symbirth=100 * s.own_asym.mean(), mwu_p_high_vs_low=mw.pvalue if mw else np.nan))
pd.DataFrame(rows).to_csv('lyso_daughters_summary.csv', index=False)
print('done')
