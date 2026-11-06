import re
import numpy as np
import pandas as pd
from .core import thin

KAPPA = re.compile(r"^(IGKC|IGKV)")
LAMBDA = re.compile(r"^IGL[CV]")
VGENE = re.compile(r"^IG[HKL]V")


def _stats(M, genes):
    k = np.array([bool(KAPPA.match(g)) for g in genes])
    l = np.array([bool(LAMBDA.match(g)) for g in genes])
    v = np.array([bool(VGENE.match(g)) for g in genes])
    hk = np.asarray(M[:, k].sum(axis=1)).ravel() > 0 if k.any() else np.zeros(M.shape[0], bool)
    hl = np.asarray(M[:, l].sum(axis=1)).ravel() > 0 if l.any() else np.zeros(M.shape[0], bool)
    nv = np.asarray((M[:, v] > 0).sum(axis=1)).ravel() if v.any() else np.zeros(M.shape[0])
    ig = hk | hl
    return dict(n=M.shape[0], pct_Ig_positive=round(100 * ig.mean(), 1),
                pct_kappa_and_lambda_of_Ig=round(100 * (hk & hl).sum() / max(ig.sum(), 1), 1),
                median_V_genes_in_Ig_cells=float(np.median(nv[ig])) if ig.any() else 0.0)


def ig_clonality(dm, groups, plasma_type=None):
    """A real B or plasma cell expresses one light chain and one V gene per chain; ambient antibody RNA
    pooled from many clones shows kappa with lambda and many V genes in the same cell."""
    rows = []
    for g in groups:
        qr = dm.query_rows(g)
        rows.append(dict(set=g, **_stats(dm.C[qr], dm.genes)))
        if plasma_type:
            pr = dm._samp(dm.pool(plasma_type), dm.n_ds)
            rows.append(dict(set=f"real {plasma_type} at depth of {g}", **_stats(thin(dm.C[pr], dm.rng.choice(dm.depth[qr], len(pr)), dm.rng), dm.genes)))
    return pd.DataFrame(rows).set_index("set")
