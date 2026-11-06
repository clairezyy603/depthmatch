import argparse
import os
import yaml
import scanpy as sc
from .core import DepthMatch
from .mixed import split_mixed
from .clonality import ig_clonality


def main(argv=None):
    p = argparse.ArgumentParser(prog="depthmatch", description="Depth-matched identity calls for low-information single-cell groups")
    p.add_argument("config", help="YAML config, see examples/config_example.yaml")
    a = p.parse_args(argv)
    cfg = yaml.safe_load(open(a.config))
    out = cfg.get("out", "depthmatch_results")
    adata = sc.read_h5ad(cfg["h5ad"])
    dm = DepthMatch(adata, cfg["label_key"], cfg["query_key"], layer=cfg.get("layer", "counts"), batch_key=cfg.get("batch_key"),
                    seed=cfg.get("seed", 0), **cfg.get("params", {}))
    res = dm.run(cfg["candidates"], cfg.get("reference_types"))
    res.save(out)
    print(res.summary.to_string())
    for m in cfg.get("mixed", []):
        per_cell, summ = split_mixed(dm, m["group"], m["type_a"], m["type_b"])
        tag = f"mixed_{m['group']}_{m['type_a']}_{m['type_b']}".replace(" ", "_").replace("/", "_")
        per_cell.to_csv(os.path.join(out, f"{tag}_per_cell.csv"))
        summ.to_csv(os.path.join(out, f"{tag}_summary.csv"))
        print(f"\n{m['group']}: {m['type_a']} vs {m['type_b']}\n" + summ.to_string())
    if cfg.get("clonality"):
        c = cfg["clonality"]
        t = ig_clonality(dm, c["groups"], c.get("plasma_type"))
        t.to_csv(os.path.join(out, "ig_clonality.csv"))
        print("\n" + t.to_string())
    print(f"\nwrote results to {out}/")


if __name__ == "__main__":
    main()
