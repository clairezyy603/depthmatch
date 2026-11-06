import numpy as np
import pandas as pd
from .core import thin


def _frac(M, idx):
    tot = np.asarray(M.sum(axis=1)).ravel()
    return np.asarray(M[:, idx].sum(axis=1)).ravel() / np.maximum(tot, 1)


def split_mixed(dm, group, type_a, type_b, n_doublets=1500, percentile=95, n_program=30, min_own=0.2, max_other=0.1):
    """Per-cell test of whether a group carrying two programs holds doublets or a faint second signal.

    Builds a program for each type from full-depth reference cells. A program counts as present above the given
    percentile of the opposite pure type. Cells with both programs are called doublet-like when both scores also reach
    the lower tail of synthetic doublets of the two types, and faint second program otherwise.
    Returns (per_cell, summary).
    """
    rng = dm.rng
    qr = dm.query_rows(group)
    Mq, depth = dm.C[qr], dm.depth[qr]
    ra, rb = dm._samp(dm.pool(type_a), dm.n_ref), dm._samp(dm.pool(type_b), dm.n_ref)
    up_a, up_b = dm._de(dm._full(ra), dm._full(rb)), dm._de(dm._full(rb), dm._full(ra))
    Ra, Rb = thin(dm.C[ra], rng.choice(depth, len(ra)), rng), thin(dm.C[rb], rng.choice(depth, len(rb)), rng)
    gi = {g: i for i, g in enumerate(dm.genes)}

    def program(gl, own, other):
        idx = [gi[g] for g in gl]
        if not idx:
            return []
        do = np.asarray((own[:, idx] > 0).mean(axis=0)).ravel()
        dt = np.asarray((other[:, idx] > 0).mean(axis=0)).ravel()
        return [g for g, a, b in zip(gl, do, dt) if a >= min_own and b <= max_other][:n_program]

    pa, pb = program(up_a, Ra, Rb), program(up_b, Rb, Ra)
    if not pa or not pb:
        raise ValueError(f"could not build both programs ({type_a}: {len(pa)} genes, {type_b}: {len(pb)} genes)")
    ia, ib = [gi[g] for g in pa], [gi[g] for g in pb]
    da, db = rng.choice(len(ra), n_doublets), rng.choice(len(rb), n_doublets)
    D = thin(dm.C[ra[da]] + dm.C[rb[db]], rng.choice(depth, n_doublets), rng)
    sets = {group: Mq, f"real {type_a}": Ra, f"real {type_b}": Rb, f"synthetic {type_a}+{type_b} doublet": D}
    F = {n: (_frac(M, ia), _frac(M, ib)) for n, M in sets.items()}
    ta = np.percentile(F[f"real {type_b}"][0], percentile)
    tb = np.percentile(F[f"real {type_a}"][1], percentile)

    dn = f"synthetic {type_a}+{type_b} doublet"
    la, lb = np.percentile(F[dn][0], 100 - percentile), np.percentile(F[dn][1], 100 - percentile)
    classes = ["doublet-like", "faint second program", f"{type_a} only", f"{type_b} only", "neither"]

    def classify(fa, fb):
        a, b = fa > ta, fb > tb
        dbl = (fa >= la) & (fb >= lb)
        both = np.where(dbl, classes[0], classes[1])
        return np.where(a & b, both, np.where(a, classes[2], np.where(b, classes[3], classes[4])))

    rows = []
    for n, (fa, fb) in F.items():
        v = pd.Series(classify(fa, fb)).value_counts(normalize=True).mul(100).round(1)
        rows.append(dict(set=n, n=len(fa), **{c: v.get(c, 0.0) for c in classes},
                         **{f"median_score_{type_a}": round(float(np.median(fa)), 4), f"median_score_{type_b}": round(float(np.median(fb)), 4)}))
    summary = pd.DataFrame(rows).set_index("set")
    summary.attrs.update(program_a=pa, program_b=pb, threshold_a=ta, threshold_b=tb, doublet_floor_a=la, doublet_floor_b=lb)
    fa, fb = F[group]
    per_cell = pd.DataFrame({f"score_{type_a}": fa, f"score_{type_b}": fb, "call": classify(fa, fb)}, index=dm.obs_names[qr])
    return per_cell, summary
