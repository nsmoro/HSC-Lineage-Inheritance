"""Step 3 - analyses requested after the first presentation.

Usage (inside the results folder, after step 1):
  python 03_lineage_analyses.py basic           -> n_numbers.csv, pctdiff_by_birth.csv (+ pctdiff_raw.pkl),
                                                   generation_trends.csv
  python 03_lineage_analyses.py patterns_mouse  -> patterns_mouse.csv  (lineages of 3 divisions)
  python 03_lineage_analyses.py patterns_human  -> patterns_human.csv  (lineages of 2 divisions)
  python 03_lineage_analyses.py sistcous        -> sisters_cousins.csv

  n_numbers          experiments, trees, cells and divisions per generation
  pctdiff_by_birth   threshold-free: % difference of a cell's own division by how it was born; shift in
                     percentile points (born asymmetric - born symmetric) with 95% CI (colony resampling)
  generation_trends  median marker level, absolute and % difference per daughter generation
  patterns           every lineage written as one letter per division (A = more unequal than the median
                     division of its stratum, S = more equal); observed counts vs permutation expectation
  sisters_cousins    % of sister / cousin pairs dividing in the same class (thirds), by founder division class
"""
import sys, os, itertools, warnings; warnings.filterwarnings('ignore'); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
from dataset import KEY
rng = np.random.default_rng(11)
MEAS = ['log2ratio', 'pctdiff', 'absdiff', 'absdiff_adj']
cells, div, pairs = pd.read_pickle('dataset.pkl')
STRATA = ['Experiment', 'Condition', 'Marker', 'DaughterGen']

def groups_by_colony(df):
    ids = (df.Experiment.astype(str) + '|' + df.Condition.astype(str) + '|' + df.Position.astype(str) + '|' + df.Colony.astype(str)).values
    u, inv = np.unique(ids, return_inverse=True); return [np.flatnonzero(inv == i) for i in range(len(u))]

def n_numbers():
    rows = []
    for sp, c in cells.groupby('Species'):
        c = c.copy(); c['tree'] = c.Experiment + '|' + c.Condition + '|' + c.Position.astype(str) + '|' + c.Colony.astype(str)
        dv = c[c.Fate == 'divided'].groupby('tree').size()
        r = dict(Species=sp, experiments=c.Experiment.nunique(), conditions=c.groupby('Experiment').Condition.nunique().sum(),
                 positions=c.groupby(['Experiment', 'Condition']).Position.nunique().sum() if sp == 'human' else c.groupby('Experiment').Position.nunique().sum(),
                 trees_tracked=c.tree.nunique(), trees_with_division=int((dv > 0).sum()), cells=len(c))
        for g in range(0, 5):
            r[f'cells_gen{g}'] = int((c.Generation == g).sum()); r[f'divided_gen{g}'] = int(((c.Generation == g) & (c.Fate == 'divided')).sum())
        r['max_depth_tree_gen'] = int(c.Generation.max())
        rows.append(r)
    return pd.DataFrame(rows)

def pctdiff_by_birth():
    """Continuous: % difference of a cell's own division, by how it was born (parent class, tertiles of pctdiff)."""
    p = pairs.dropna(subset=['p_pctdiff_class', 'c_pctdiff', 'c_pctdiff_pct']).copy()
    p['born'] = p.p_pctdiff_class.astype(str)
    p.loc[p.born == 'asym', 'born'] = np.where(p.loc[p.born == 'asym', 'mother_was_high'], 'asym_more', 'asym_less')
    rows = []
    for (sp, mk), g in p.groupby(['c_Species', 'Marker']):
        grp = groups_by_colony(g); n = len(grp)
        def shift(x):  # mean percentile (within stratum) born asym minus born sym
            return x[x.born.str.startswith('asym')].c_pctdiff_pct.mean() - x[x.born == 'sym'].c_pctdiff_pct.mean()
        bs = [shift(g.iloc[np.concatenate([grp[i] for i in rng.integers(0, n, n)])]) for _ in range(200)]
        r = dict(Species=sp, Marker=mk, shift_pct_points=100 * shift(g), lo=100 * np.nanpercentile(bs, 2.5), hi=100 * np.nanpercentile(bs, 97.5))
        for b in ['sym', 'int', 'asym_more', 'asym_less']:
            v = g[g.born == b].c_pctdiff
            r[f'median_{b}'] = v.median(); r[f'n_{b}'] = len(v)
        rows.append(r)
    return pd.DataFrame(rows), p[['c_Species', 'Marker', 'born', 'c_pctdiff', 'c_DaughterGen']]

def generation_trends():
    d = div.copy(); d['total'] = d.v_high + d.v_low
    rows = []
    for (sp, mk, g), h in d.groupby(['Species', 'Marker', 'DaughterGen']):
        if len(h) < 20: continue
        rows.append(dict(Species=sp, Marker=mk, DaughterGen=g, n=len(h), absdiff=h.absdiff.median(), total=h.total.median(),
                         pctdiff=h.pctdiff.median(), absdiff_q1=h.absdiff.quantile(.25), absdiff_q3=h.absdiff.quantile(.75),
                         pctdiff_q1=h.pctdiff.quantile(.25), pctdiff_q3=h.pctdiff.quantile(.75)))
    return pd.DataFrame(rows)

def pattern_counts(length, species, n_perm=200):
    """Lineages of `length` consecutive divisions; each division A (more unequal than the median of its stratum) or
    S (more equal). Count each pattern; expected = permutation within experiment x condition x marker x generation."""
    out = []
    for me in MEAS:
        d = div[(div.Species == species)].dropna(subset=[f'{me}_pct']).copy()
        d['gen'] = np.floor(np.log2(d.Mother)).astype(int); d = d[d.gen < length]
        d['A'] = (d[f'{me}_pct'] > 0.5).values
        lut = {k: i for i, k in enumerate(zip(d.Experiment, d.Condition, d.Position, d.Colony, d.Marker, d.Mother))}
        last = d[d.gen == length - 1]; idx, mks = [], []
        for r in last.itertuples(index=False):
            ks = [(r.Experiment, r.Condition, r.Position, r.Colony, r.Marker, r.Mother >> s) for s in range(length - 1, -1, -1)]
            if all(k in lut for k in ks): idx.append([lut[k] for k in ks]); mks.append(r.Marker)
        if not idx: continue
        idx = np.array(idx); mks = np.array(mks)
        strata = d.groupby(['Experiment', 'Condition', 'Marker', 'gen']).ngroup().values
        grs = [np.flatnonzero(strata == s) for s in np.unique(strata)]
        A = d.A.values; pats = [''.join(p) for p in itertools.product('AS', repeat=length)]
        code = lambda Av: np.array([''.join('A' if v else 'S' for v in row) for row in Av[idx]])
        obs = code(A); nulls = []
        for _ in range(n_perm):
            Ap = A.copy()
            for gr in grs: Ap[gr] = A[rng.permutation(gr)]
            nulls.append(code(Ap))
        for mk in np.unique(mks):
            sel = mks == mk; n = sel.sum()
            for p in pats:
                o = (obs[sel] == p).sum(); e = np.array([(nl[sel] == p).sum() for nl in nulls])
                pv = (1 + np.sum(np.abs(e - e.mean()) >= abs(o - e.mean()))) / (n_perm + 1)
                out.append(dict(Species=species, Measure=me, Marker=mk, length=length, pattern=p, n_lineages=int(n), count=int(o),
                                pct=100 * o / n, expected_count=e.mean(), expected_pct=100 * e.mean() / n, p=pv))
    return pd.DataFrame(out)

def sisters_cousins(n_boot=200):
    """Similarity of sisters' and cousins' own divisions, split by the founder division (cell 1) class.
    Similarity = % of pairs in the same class (S/I/A tertiles; chance 33%)."""
    rows = []
    for me in MEAS:
        d = div.dropna(subset=[f'{me}_class']).copy()
        key = ['Experiment', 'Condition', 'Position', 'Colony', 'Marker']
        cl = d.set_index(key + ['Mother'])[f'{me}_class']
        founder = d[d.Mother == 1].set_index(key)[f'{me}_class'].rename('founder')
        pairs_def = {'sisters (cells 2 vs 3)': [(2, 3)], 'sisters (gen 2: 4v5, 6v7)': [(4, 5), (6, 7)], 'cousins (4/5 vs 6/7)': [(4, 6), (4, 7), (5, 6), (5, 7)]}
        for name, plist in pairs_def.items():
            recs = []
            for a, b in plist:
                A = d[d.Mother == a].set_index(key)[f'{me}_class']; B = d[d.Mother == b].set_index(key)[f'{me}_class']
                j = pd.concat([A.rename('a'), B.rename('b')], axis=1, join='inner').join(founder, how='inner')
                recs.append(j)
            j = pd.concat(recs).reset_index()
            if j.empty: continue
            j['same'] = (j.a.astype(str) == j.b.astype(str)).astype(float)
            sp_of = cells[['Experiment', 'Species']].drop_duplicates().set_index('Experiment').Species
            j['Species'] = j.Experiment.map(sp_of)
            for (sp, mk, f), g in j.groupby(['Species', 'Marker', 'founder']):
                if len(g) < 15: continue
                grp = groups_by_colony(g); n = len(grp)
                bs = [g.iloc[np.concatenate([grp[i] for i in rng.integers(0, n, n)])].same.mean() for _ in range(n_boot)]
                rows.append(dict(Species=sp, Measure=me, Marker=mk, relation=name, founder=f, n_pairs=len(g),
                                 pct_same=100 * g.same.mean(), lo=100 * np.nanpercentile(bs, 2.5), hi=100 * np.nanpercentile(bs, 97.5)))
    return pd.DataFrame(rows)

if __name__ == '__main__':
    w = sys.argv[1]
    if w == 'basic':
        n_numbers().to_csv('n_numbers.csv', index=False)
        a, raw = pctdiff_by_birth(); a.to_csv('pctdiff_by_birth.csv', index=False); raw.to_pickle('pctdiff_raw.pkl')
        generation_trends().to_csv('generation_trends.csv', index=False)
    elif w == 'patterns_mouse':
        pattern_counts(3, 'mouse').to_csv('patterns_mouse.csv', index=False)
    elif w == 'patterns_human':
        pattern_counts(2, 'human').to_csv('patterns_human.csv', index=False)
    elif w == 'sistcous':
        sisters_cousins().to_csv('sisters_cousins.csv', index=False)
    print('done', w)
