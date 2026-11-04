"""Step 2 - do daughters pass on what they received?

For every division (each species, measure and marker) we ask how often a cell divides asymmetrically itself,
depending on how it was born: in a symmetric, intermediate or asymmetric division (thirds of the measure within
experiment x condition x marker x generation), and - if born asymmetric - whether it received more or less.

Usage (inside the results folder, after step 1):
  python 02_birth_vs_division.py received   -> birth_vs_division.csv
        % dividing asymmetrically per birth group, 95% CI by resampling colonies, and the odds ratio
        (born asymmetric vs born symmetric) with its CI
  python 02_birth_vs_division.py lin3       -> lineage_enrichment_3div_mouse.csv
  python 02_birth_vs_division.py lin2       -> lineage_enrichment_2div.csv
        lineages of 3 (mouse) or 2 consecutive divisions: observed vs random frequency of lineages that are
        asymmetric in every division, never asymmetric, and the number of switches (permutation of division
        classes within experiment x condition x marker x generation, 200 permutations)
"""
import sys, os, warnings; warnings.filterwarnings('ignore'); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd, itertools
from dataset import KEY
rng = np.random.default_rng(7)
MEAS = ['log2ratio', 'pctdiff', 'absdiff', 'absdiff_adj']

def col_groups(df):
    ids = (df.Experiment.astype(str) + '|' + df.Condition.astype(str) + '|' + df.Position.astype(str) + '|' + df.Colony.astype(str)).values
    u, inv = np.unique(ids, return_inverse=True)
    return [np.flatnonzero(inv == i) for i in range(len(u))]

def received_vs_division(pairs, n_boot=200):
    """For each species/measure/marker: % of cells that divide asymmetrically, by how they were born:
    symmetric, intermediate, asymmetric-received-more, asymmetric-received-less."""
    out = []
    for (sp, mk), g in pairs.groupby(['c_Species', 'Marker']):
        for me in MEAS:
            h = g.dropna(subset=[f'p_{me}_class', f'c_{me}_class']).copy()
            if len(h) < 30: continue
            h['born'] = h[f'p_{me}_class'].astype(str)
            h.loc[h.born == 'asym', 'born'] = np.where(h.loc[h.born == 'asym', 'mother_was_high'], 'asym_more', 'asym_less')
            h['div_asym'] = (h[f'c_{me}_class'] == 'asym').astype(float)
            grp = col_groups(h); n = len(grp)
            boots = []
            for _ in range(n_boot):
                idx = np.concatenate([grp[i] for i in rng.integers(0, n, n)]); b = h.iloc[idx]
                r = b.groupby('born').div_asym.mean()
                A = b[b.born.str.startswith('asym')].div_asym; S = b[b.born == 'sym'].div_asym
                pa, ps = (A.sum() + .5) / (len(A) + 1), (S.sum() + .5) / (len(S) + 1)
                boots.append({**r.to_dict(), 'log2OR': np.log2(pa / (1 - pa) / (ps / (1 - ps)))})
            B = pd.DataFrame(boots)
            obs = h.groupby('born').div_asym.mean()
            A = h[h.born.str.startswith('asym')].div_asym; S = h[h.born == 'sym'].div_asym
            pa, ps = (A.sum() + .5) / (len(A) + 1), (S.sum() + .5) / (len(S) + 1)
            obs['log2OR'] = np.log2(pa / (1 - pa) / (ps / (1 - ps)))
            for k in ['sym', 'int', 'asym_more', 'asym_less', 'log2OR']:
                if k not in obs or (k != 'log2OR' and (h.born == k).sum() < 15): continue
                out.append(dict(Species=sp, Measure=me, Marker=mk, group=k, n=int((h.born == k).sum()) if k != 'log2OR' else len(h),
                                value=obs[k], lo=np.nanpercentile(B[k], 2.5), hi=np.nanpercentile(B[k], 97.5)))
    return pd.DataFrame(out)

def lineage_enrichment(div, chains_len, n_perm=200):
    """Lineages of `chains_len` consecutive divisions: observed/expected (random) for uniformly asymmetric (all A),
    never asymmetric (all N) and number of switches, per species/measure/marker."""
    res = []
    for me in MEAS:
        d = div.dropna(subset=[f'{me}_class']).copy()
        d['gen'] = np.floor(np.log2(d.Mother)).astype(int); d = d[d.gen < chains_len]
        d['A'] = (d[f'{me}_class'] == 'asym').values
        lut = {k: i for i, k in enumerate(zip(d.Experiment, d.Condition, d.Position, d.Colony, d.Marker, d.Mother))}
        last = d[d.gen == chains_len - 1]
        idx, meta = [], []
        for r in last.itertuples(index=False):
            ks = [(r.Experiment, r.Condition, r.Position, r.Colony, r.Marker, r.Mother >> s) for s in range(chains_len - 1, -1, -1)]
            if all(k in lut for k in ks):
                idx.append([lut[k] for k in ks]); meta.append((r.Species, r.Marker))
        if not idx: continue
        idx = np.array(idx); meta = pd.DataFrame(meta, columns=['Species', 'Marker'])
        strata = d.groupby(['Experiment', 'Condition', 'Marker', 'gen']).ngroup().values
        groups = [np.flatnonzero(strata == s) for s in np.unique(strata)]
        A = d.A.values
        def summ(Av):
            p = Av[idx]; return p.all(1), (~p).all(1), (p[:, 1:] != p[:, :-1]).sum(1)
        oA, oN, oS = summ(A)
        nulls = []
        for _ in range(n_perm):
            Ap = A.copy()
            for gr in groups: Ap[gr] = A[rng.permutation(gr)]
            nulls.append(summ(Ap))
        for (sp, mk), g in meta.groupby(['Species', 'Marker']):
            rows = g.index.values
            for name, o, k in [('all_asym', oA, 0), ('never_asym', oN, 1), ('switches', oS, 2)]:
                ob = o[rows].mean(); nl = np.array([n[k][rows].mean() for n in nulls])
                p = (1 + np.sum(np.abs(nl - nl.mean()) >= abs(ob - nl.mean()))) / (n_perm + 1)
                res.append(dict(Species=sp, Measure=me, Marker=mk, chain_len=chains_len, n_lineages=len(rows), stat=name,
                                observed=ob, expected=nl.mean(), ratio=ob / nl.mean() if nl.mean() > 0 else np.nan, p=p))
    return pd.DataFrame(res)

if __name__ == '__main__':
    cells, div, pairs = pd.read_pickle('dataset.pkl')
    which = sys.argv[1]
    if which == 'received':
        received_vs_division(pairs).to_csv('birth_vs_division.csv', index=False)
    elif which == 'lin3':
        lineage_enrichment(div[div.Species == 'mouse'], 3).to_csv('lineage_enrichment_3div_mouse.csv', index=False)
    elif which == 'lin2':
        lineage_enrichment(div, 2).to_csv('lineage_enrichment_2div.csv', index=False)
    print('done', which)
