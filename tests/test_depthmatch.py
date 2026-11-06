import depthmatch as dm
from depthmatch.simulate import simulate


def setup():
    a = simulate(n_per_type=400, n_low=150, seed=1)
    return dm.DepthMatch(a, "cell_type", "query_group", batch_key="batch", n_ref=400, n_ds=300, n_cal=150)


def test_identity_calls():
    d = setup()
    res = d.run({"lowMac": ["Macrophage", "Plasma"], "lowAT2": ["AT2", "Macrophage"]})
    assert res.summary.loc["lowMac", "call"] == "Macrophage"
    assert res.summary.loc["lowAT2", "call"] == "AT2"


def test_mixed_signal():
    d = setup()
    _, dbl = dm.split_mixed(d, "mixMacAT2", "AT2", "Macrophage", n_doublets=500)
    _, amb = dm.split_mixed(d, "ambAT2", "AT2", "Macrophage", n_doublets=500)
    assert dbl.loc["mixMacAT2", "doublet-like"] > 60
    assert amb.loc["ambAT2", "doublet-like"] < 10


def test_clonality():
    d = setup()
    t = dm.ig_clonality(d, ["lowMac"], plasma_type="Plasma")
    assert t.loc["lowMac", "pct_kappa_and_lambda_of_Ig"] > 50
    assert t.iloc[1]["pct_kappa_and_lambda_of_Ig"] < 10
