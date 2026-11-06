# depthmatch

Depth-matched identity calls for low-information single-cell clusters.

Most single-cell atlases end up with a cluster of cells that have very few genes. Ambient RNA tools (CellBender, SoupX, DecontX) remove contamination and doublet tools flag doublets, but neither answers the question an annotator actually faces: **what are these cells?** Comparing them directly with normal cells does not work, because at low depth almost every gene looks "missing", so any reference type looks different.

`depthmatch` makes that comparison fair. It downsamples real reference cells to the same sequencing depth as the low-information cells, measures how different two groups of the *same* type look at that depth, and only then asks which reference type the query matches.

## How it works

1. **Identity genes.** For every reference type, genes that are up against every other type are found with pairwise Wilcoxon tests on full-depth, high-quality reference cells (cells at or above their type's median gene count). Only genes that are still detected in at least 30% of that type's cells after downsampling to the query depth are kept, so the test never relies on genes that cannot be seen at low depth.
2. **Calibration.** Each reference type is split in two (by batch when a batch key is given), both halves are downsampled to the query depth, and the number of differentially expressed genes between them is recorded. This is the noise floor: how many genes two identical populations differ by at this depth. The same is done between different types.
3. **Comparison.** Each query group is compared in both directions with each candidate type downsampled to that group's depth distribution. A candidate **matches** when none of its identity genes are up in the reference relative to the query, is **near** when one or two are, and is **not** a match otherwise. The calibrated call reports whether the number of differing genes stays within the same-type noise floor.

Two optional checks handle the usual complications:

- **Mixed signals** (`split_mixed`). For a group that expresses two programs, each cell is scored for both programs and compared with pure reference cells and with **synthetic doublets** made from real cells of the two types. Cells are called `doublet-like` when both programs reach doublet levels, and `faint second program` when the second program is present but far weaker, which is the pattern ambient RNA produces.
- **Antibody clonality** (`ig_clonality`). A real B or plasma cell expresses one light chain and one V gene per chain. Ambient antibody RNA pooled from many clones shows kappa together with lambda and many V genes in the same cell. The check compares each group with real plasma cells downsampled to the same depth.

## Install

```bash
git clone https://github.com/clairezyy603/depthmatch.git
cd depthmatch
conda env create -f environment.yml
conda activate depthmatch
```

or, in an existing environment with scanpy: `pip install -e .`

## Quick start

Your AnnData needs raw counts (in a layer or in `X`), a column of reference cell-type labels, and a column that marks the low-information groups to test (empty or NaN for reference cells).

```python
import scanpy as sc
import depthmatch as dm

adata = sc.read_h5ad("my_data.h5ad")
d = dm.DepthMatch(adata, label_key="cell_type", query_key="low_quality_group",
                  layer="counts", batch_key="library_id")

res = d.run({"groupA": ["Alveolar macrophage", "Monocyte", "Plasma cell"],
             "groupB": ["AT2", "Alveolar macrophage"]})
print(res.summary)
res.save("depthmatch_results")

per_cell, summary = dm.split_mixed(d, "groupB", "AT2", "Alveolar macrophage")
print(dm.ig_clonality(d, ["groupA"], plasma_type="Plasma cell"))
```

From the command line, with a config file (see `examples/config_example.yaml`):

```bash
depthmatch examples/config_example.yaml
```

## Outputs

| File | Contents |
|---|---|
| `summary.csv` | One call per query group, with the result for every candidate |
| `comparisons.csv` | Every group by candidate comparison: genes up in each direction, identity genes missing, calibrated call |
| `identity_genes.csv` | Identity genes used for each reference type and the rule that selected them |
| `calibration_same_type.csv` | Genes separating two halves of the same type at query depth (the noise floor) |
| `calibration_between_types.csv` | Genes separating each pair of types at query depth |
| `mixed_*_summary.csv` | Class percentages and median program scores for the group, both pure references and synthetic doublets |
| `ig_clonality.csv` | Kappa and lambda co-expression and V genes per cell, for each group and for real plasma cells at the same depth |

**Reading the mixed-signal table.** Read the query row against the other three. The synthetic doublet row shows how often real doublets are detected at this depth, and the pure reference rows show the false-positive rate. A group of doublets looks like the synthetic doublet row, while a group with ambient RNA has mostly `faint second program` cells and a median second-program score far below the doublets.

## Demo

`examples/simulated_demo.py` builds a toy dataset with four reference types and four low-information groups whose truth is known (shallow macrophages carrying ambient antibody RNA, shallow AT2 cells, true doublets, and AT2 cells with a faint macrophage signal) and runs every module. It finishes in about ten seconds:

```
group   call
lowMac  Macrophage   (Plasma: not, T cell: not)
lowAT2  AT2          (Macrophage: not)
mixMacAT2  83% doublet-like   (synthetic doublets 90%)
ambAT2      0% doublet-like, 70% faint second program
lowMac  kappa + lambda in 95% of Ig+ cells   (real plasma 0.2%)
```

Run the tests with `pytest tests`.

## Limits

- The reference labels must already be trustworthy, and the right type has to be among the candidates; the tool tests hypotheses, it does not discover new types.
- Identity genes need at least two reference types, and the identity call needs at least three detectable identity genes; otherwise it reports `too few identity genes`.
- The thresholds (30% detection, 0 / 1–2 / 3+ missing identity genes, the 95th percentile for programs) are conventions. Read the gene lists and counts, not only the call.
- Downsampling matches total UMIs, not gene counts, so very different library complexity between datasets can still leave some depth difference.

## Origin

Developed while annotating a cluster of about 36,000 low-gene-count cells in a multi-dataset lung single-cell atlas, where it resolved the cluster into macrophages carrying ambient RNA, AT2 cells and a small doublet group.

## License

MIT
