import os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import scipy.sparse as sp
from collections import defaultdict

EXCLUDE = r"^(RPL|RPS|MT-|MTRNR2L)"


def counts_matrix(adata, layer=None):
    M = adata.layers[layer] if layer else adata.X
    M = sp.csr_matrix(M, dtype=np.float32)
    M.data = np.rint(M.data)
    return M


def thin(M, targets, rng):
    """Binomially downsample each row of a count matrix to a target UMI total."""
    M = sp.csr_matrix(M, dtype=np.float32, copy=True)
    s = np.asarray(M.sum(axis=1)).ravel()
    p = np.minimum(1.0, np.asarray(targets, dtype=float) / np.maximum(s, 1))
    r = np.repeat(np.arange(M.shape[0]), np.diff(M.indptr))
    M.data = rng.binomial(M.data.astype(np.int64), p[r]).astype(np.float32)
    M.eliminate_zeros()
    return M


def to_adata(M, genes):
    x = ad.AnnData(sp.csr_matrix(M, dtype=np.float32), var=pd.DataFrame(index=genes))
    x.obs_names = [f"c{i}" for i in range(x.n_obs)]
    x.layers["counts"] = x.X.copy()
    sc.pp.normalize_total(x, target_sum=1e4)
    sc.pp.log1p(x)
    return x


def markers(X, Y, min_lfc=1.0, min_in=0.25, max_out=0.25, max_padj=0.05, n=200, exclude=EXCLUDE):
    """Genes up in X versus Y, tested on the two groups only (Wilcoxon)."""
    fx = np.asarray((X.layers["counts"] > 0).mean(axis=0)).ravel()
    fy = np.asarray((Y.layers["counts"] > 0).mean(axis=0)).ravel()
    keep = (fx >= min_in) & (fy <= max_out)
    if exclude:
        keep &= ~X.var_names.str.match(exclude)
    if not keep.any():
        return []
    g = X.var_names[keep]
    z = ad.concat([X[:, g], Y[:, g]], index_unique="-")
    z.obs["grp"] = pd.Categorical(["x"] * X.n_obs + ["y"] * Y.n_obs, categories=["x", "y"])
    sc.tl.rank_genes_groups(z, "grp", groups=["x"], reference="y", method="wilcoxon")
    d = sc.get.rank_genes_groups_df(z, "x")
    d = d[(d.logfoldchanges >= min_lfc) & (d.pvals_adj <= max_padj)]
    return list(d.sort_values("scores", ascending=False)["names"].astype(str)[:n])


class Result:
    def __init__(self, **tables):
        self.tables = tables
        for k, v in tables.items():
            setattr(self, k, v)

    def save(self, out):
        os.makedirs(out, exist_ok=True)
        for k, v in self.tables.items():
            v.to_csv(os.path.join(out, f"{k}.csv"))


class DepthMatch:
    """Depth-matched identity calls for low-information cell groups.

    adata: counts in `layer` (or X); `label_key`: reference cell-type labels;
    `query_key`: the low-information groups to test (empty or NaN for reference cells).
    """

    def __init__(self, adata, label_key, query_key, layer="counts", batch_key=None, seed=0,
                 n_ref=1500, n_ds=600, n_cal=300, n_query=3000, n_splits=3, min_detect=0.3, n_identity=20,
                 min_lfc=1.0, min_in=0.25, max_out=0.25, max_padj=0.05, exclude=EXCLUDE):
        self.rng = np.random.default_rng(seed)
        self.obs_names = np.asarray(adata.obs_names)
        self.C = counts_matrix(adata, layer)
        self.genes = adata.var_names.copy()
        self.labels = adata.obs[label_key].astype(str).values
        q = adata.obs[query_key]
        self.query = np.where(q.isna(), "", q.astype(str)).astype(object)
        self.query[np.isin(self.query, ["nan", "None"])] = ""
        self.batch = adata.obs[batch_key].astype(str).values if batch_key else None
        self.depth = np.asarray(self.C.sum(axis=1)).ravel()
        self.n_genes = np.asarray((self.C > 0).sum(axis=1)).ravel()
        self.n_ref, self.n_ds, self.n_cal, self.n_query, self.n_splits = n_ref, n_ds, n_cal, n_query, n_splits
        self.min_detect, self.n_identity = min_detect, n_identity
        self.de_kw = dict(min_lfc=min_lfc, min_in=min_in, max_out=max_out, max_padj=max_padj, exclude=exclude)
        self.pools, self.identity = {}, {}

    def _samp(self, idx, n):
        return idx if len(idx) <= n else np.sort(self.rng.choice(idx, n, replace=False))

    def _full(self, rows):
        return to_adata(self.C[rows], self.genes)

    def _ds(self, rows, target_depths):
        return to_adata(thin(self.C[rows], self.rng.choice(target_depths, len(rows)), self.rng), self.genes)

    def _de(self, X, Y):
        return markers(X, Y, **self.de_kw)

    def query_rows(self, group):
        return np.where(self.query == group)[0]

    def pool(self, t):
        if t not in self.pools:
            idx = np.where((self.labels == t) & (self.query == ""))[0]
            if len(idx) == 0:
                raise ValueError(f"no reference cells labelled {t!r}")
            self.pools[t] = idx[self.n_genes[idx] >= np.median(self.n_genes[idx])]
        return self.pools[t]

    def identity_genes(self, types, target_depths):
        hq = {t: self._full(self._samp(self.pool(t), self.n_ref)) for t in types}
        self._dsc = {t: self._ds(self._samp(self.pool(t), self.n_ds), target_depths) for t in types}
        cnt, rk = {t: defaultdict(int) for t in types}, {t: defaultdict(list) for t in types}
        for x in types:
            for y in types:
                if x == y:
                    continue
                for i, g in enumerate(self._de(hq[x], hq[y])):
                    cnt[x][g] += 1
                    rk[x][g].append(i)
        rows = []
        for x in types:
            need = len(types) - 1
            c = pd.Series(cnt[x], dtype=float)
            if c.empty:
                self.identity[x] = []
                rows.append(dict(type=x, rule="none", n_candidates=0, n_detectable=0, genes=""))
                continue
            r = pd.Series({g: np.mean(v) for g, v in rk[x].items()})
            gl, rule = c[c == need].index, "up vs all other types"
            if len(gl) < 3:
                gl, rule = c[c >= need / 2].index, "up vs at least half of other types"
            det = pd.Series(np.asarray((self._dsc[x][:, list(gl)].layers["counts"] > 0).mean(axis=0)).ravel(), index=gl)
            o = pd.DataFrame({"n_beaten": c[gl], "mean_rank": r[gl], "detected": det}).sort_values(["n_beaten", "mean_rank"], ascending=[False, True])
            o2 = o[o.detected >= self.min_detect]
            if len(o2) < 3:
                o2 = o.sort_values("detected", ascending=False).head(max(3, len(o2)))
            self.identity[x] = list(o2.index[:self.n_identity])
            rows.append(dict(type=x, rule=rule, n_candidates=len(gl), n_detectable=int((o.detected >= self.min_detect).sum()), genes=", ".join(self.identity[x])))
        return pd.DataFrame(rows).set_index("type")

    def calibrate(self, types, target_depths):
        same = []
        for x in types:
            pool = self.pool(x)
            for s in range(self.n_splits):
                if self.batch is not None and len(set(self.batch[pool])) >= 2:
                    libs = pd.Series(self.batch[pool]).value_counts()
                    L = libs.index.values.copy()
                    self.rng.shuffle(L)
                    half = set(L[np.cumsum(libs[L].values) <= libs.sum() / 2]) or {L[0]}
                    m = np.isin(self.batch[pool], list(half))
                    how = "by batch"
                else:
                    m = self.rng.random(len(pool)) < 0.5
                    how = "random halves"
                ra, rb = pool[m], pool[~m]
                if min(len(ra), len(rb)) < 50:
                    continue
                X, Y = self._ds(self._samp(ra, self.n_cal), target_depths), self._ds(self._samp(rb, self.n_cal), target_depths)
                same.append(dict(type=x, split=s, how=how, n_ab=len(self._de(X, Y)), n_ba=len(self._de(Y, X))))
        same = pd.DataFrame(same, columns=["type", "split", "how", "n_ab", "n_ba"])
        same_max = same.assign(m=same[["n_ab", "n_ba"]].max(axis=1)).groupby("type").m.max() if len(same) else pd.Series(dtype=float)
        diff = pd.DataFrame(index=list(types), columns=list(types), dtype=float)
        for x in types:
            for y in types:
                if x != y:
                    diff.loc[x, y] = len(self._de(self._dsc[x], self._dsc[y]))
        return same, same_max, diff

    def compare(self, candidates, same_max):
        rows = []
        for g, cands in candidates.items():
            qr = self.query_rows(g)
            if len(qr) == 0:
                raise ValueError(f"no query cells in group {g!r}")
            Q = self._full(self._samp(qr, self.n_query))
            for k in cands:
                R = self._ds(self._samp(self.pool(k), self.n_ds), self.depth[qr])
                up_q, up_r = self._de(Q, R), self._de(R, Q)
                lim = same_max.get(k, np.nan)
                if np.isnan(lim):
                    cal = "no baseline"
                elif len(up_r) <= lim and len(up_q) <= lim:
                    cal = "within same-type range"
                elif len(up_r) <= lim:
                    cal = "has the program plus extra genes"
                else:
                    cal = "missing part of the program"
                idg = self.identity.get(k, [])
                miss = [x for x in idg if x in set(up_r)]
                if len(idg) < 3:
                    idc = "too few identity genes"
                elif not miss:
                    idc = "matches"
                elif len(miss) <= 2:
                    idc = "near"
                else:
                    idc = "not"
                rows.append(dict(group=g, candidate=k, n_query=len(qr), median_genes_query=int(np.median(self.n_genes[qr])),
                                 median_genes_ref_downsampled=int(np.median(np.asarray((R.layers["counts"] > 0).sum(axis=1)).ravel())),
                                 n_up_in_query=len(up_q), n_up_in_ref=len(up_r), same_type_max=lim, calibrated_call=cal,
                                 identity_genes_tested=len(idg), identity_genes_missing=", ".join(miss), identity_call=idc,
                                 up_in_query=", ".join(up_q[:15]), up_in_ref=", ".join(up_r[:15])))
        return pd.DataFrame(rows)

    @staticmethod
    def summarize(comp):
        out = []
        for g, d in comp.groupby("group", sort=False):
            d = d.assign(n_miss=d.identity_genes_missing.map(lambda s: len([x for x in s.split(", ") if x])))
            m = d[d.identity_call == "matches"]
            if len(m) == 1:
                call = m.candidate.iloc[0]
            elif len(m) > 1:
                call = "ambiguous: " + " / ".join(m.sort_values("n_up_in_ref").candidate)
            elif (d.identity_call == "near").any():
                call = "near " + d[d.identity_call == "near"].sort_values(["n_miss", "n_up_in_ref"]).candidate.iloc[0]
            else:
                call = "unresolved"
            out.append(dict(group=g, call=call, candidates_tested=", ".join(d.candidate),
                            per_candidate="; ".join(f"{r.candidate}: {r.identity_call} ({r.n_miss} missing)" for r in d.itertuples())))
        return pd.DataFrame(out).set_index("group")

    def run(self, candidates, reference_types=None):
        types = list(dict.fromkeys(reference_types or [t for v in candidates.values() for t in v]))
        qd = np.concatenate([self.depth[self.query_rows(g)] for g in candidates])
        ident = self.identity_genes(types, qd)
        same, same_max, diff = self.calibrate(types, qd)
        comp = self.compare(candidates, same_max)
        return Result(summary=self.summarize(comp), comparisons=comp.set_index(["group", "candidate"]), identity_genes=ident,
                      calibration_same_type=same.set_index(["type", "split"]), calibration_between_types=diff)
