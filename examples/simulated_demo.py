import depthmatch as dm
from depthmatch.simulate import simulate

adata = simulate()
print(adata)

d = dm.DepthMatch(adata, label_key="cell_type", query_key="query_group", layer="counts", batch_key="batch")
res = d.run({"lowMac": ["Macrophage", "Plasma", "T cell"], "lowAT2": ["AT2", "Macrophage"]},
            reference_types=["Macrophage", "Plasma", "AT2", "T cell"])
print("\nidentity calls:\n" + res.summary.to_string())
res.save("demo_results")

for g in ["mixMacAT2", "ambAT2"]:
    per_cell, summary = dm.split_mixed(d, g, "AT2", "Macrophage")
    print(f"\n{g}:\n" + summary.to_string())

print("\nantibody clonality:\n" + dm.ig_clonality(d, ["lowMac"], plasma_type="Plasma").to_string())
