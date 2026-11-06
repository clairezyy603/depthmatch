import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp
from .core import thin

IG = ["IGKC", "IGKV1-5", "IGKV3-20", "IGKV4-1", "IGLC2", "IGLV1-51", "IGLV2-14", "IGLV3-21", "IGHV3-23", "IGHV1-69", "IGHV4-34", "IGHG1"]


def simulate(n_per_type=800, n_genes=1500, n_markers=40, depth=3000, n_batches=4, n_low=300, low_depth=300, seed=0):
    """Toy data with four reference types and four low-information query groups whose truth is known:
    lowMac (shallow macrophages with ambient antibody RNA), lowAT2 (shallow AT2), mixMacAT2 (true doublets),
    ambAT2 (AT2 cells with a faint macrophage signal)."""
    rng = np.random.default_rng(seed)
    types = ["Macrophage", "Plasma", "AT2", "T cell"]
    genes = [f"G{i}" for i in range(n_genes - len(IG) - n_markers * len(types))] + IG
    genes += [f"{t.replace(' ', '')}_m{j}" for t in types for j in range(n_markers)]
    base = rng.gamma(0.6, 1.0, len(genes))
    base[[genes.index(g) for g in IG]] = 0.01
    prof = {}
    for t in types:
        w = base.copy()
        w[[genes.index(f"{t.replace(' ', '')}_m{j}") for j in range(n_markers)]] += rng.gamma(4.0, 2.0, n_markers)
        prof[t] = w / w.sum()

    def cells(t, n, d):
        P = np.tile(prof[t], (n, 1))
        if t == "Plasma":
            for i in range(n):
                light = rng.choice(["IGKC", "IGLC2"])
                v = rng.choice([g for g in IG if g.startswith("IGKV" if light == "IGKC" else "IGLV")])
                for g, w in [(light, 0.15), (v, 0.08), (rng.choice(["IGHV3-23", "IGHV1-69", "IGHV4-34"]), 0.08), ("IGHG1", 0.1)]:
                    P[i, genes.index(g)] += w
                P[i] /= P[i].sum()
        lam = P * rng.lognormal(np.log(d), 0.3, n)[:, None]
        return rng.poisson(rng.gamma(5.0, lam / 5.0))

    blocks, labels, query = [], [], []
    for t in types:
        blocks.append(cells(t, n_per_type, depth)); labels += [t] * n_per_type; query += [""] * n_per_type
    amb = np.mean(cells("Plasma", 400, depth), axis=0)
    amb = amb / amb.sum()
    low_mac = thin(sp.csr_matrix(cells("Macrophage", n_low, depth)), np.full(n_low, low_depth * 0.85), rng).toarray()
    low_mac += rng.poisson(np.outer(np.full(n_low, low_depth * 0.15), amb))
    blocks.append(low_mac); labels += ["unknown"] * n_low; query += ["lowMac"] * n_low
    blocks.append(thin(sp.csr_matrix(cells("AT2", n_low, depth)), np.full(n_low, low_depth), rng).toarray()); labels += ["unknown"] * n_low; query += ["lowAT2"] * n_low
    dbl = cells("Macrophage", n_low, depth) + cells("AT2", n_low, depth)
    blocks.append(thin(sp.csr_matrix(dbl), np.full(n_low, low_depth * 2), rng).toarray()); labels += ["unknown"] * n_low; query += ["mixMacAT2"] * n_low
    macamb = np.mean(cells("Macrophage", 400, depth), axis=0)
    at2 = thin(sp.csr_matrix(cells("AT2", n_low, depth)), np.full(n_low, low_depth * 1.9), rng).toarray()
    at2 += rng.poisson(np.outer(np.full(n_low, low_depth * 0.1), macamb / macamb.sum()))
    blocks.append(at2); labels += ["unknown"] * n_low; query += ["ambAT2"] * n_low
    X = sp.csr_matrix(np.vstack(blocks).astype(np.float32))
    obs = pd.DataFrame({"cell_type": labels, "query_group": pd.Series(query).replace("", np.nan).values,
                        "batch": [f"lib{i}" for i in rng.integers(0, n_batches, X.shape[0])]}, index=[f"cell{i}" for i in range(X.shape[0])])
    a = ad.AnnData(X, obs=obs, var=pd.DataFrame(index=genes))
    a.layers["counts"] = a.X.copy()
    return a
