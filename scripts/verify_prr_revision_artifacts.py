"""Numerical consistency checks on the saved PRR revision artifacts."""
from pathlib import Path
import json
import numpy as np
from scipy.io import loadmat

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/"data"/"prr_revision"


def csv(name):
    return np.genfromtxt(DATA/name,delimiter=",",names=True,dtype=None,encoding="utf8")


def main():
    report=json.loads((DATA/"revision_summary.json").read_text())
    grid=csv("grid_convergence.csv")
    for field in ["F_Bell","F_BIC","F_Bell_hold"]:
        assert np.ptp(grid[field])<1e-10,field
    assert grid["BIC_eigen_residual"].max()<1e-13
    assert abs(report["analyticZ"]-1/(1+0.09**2/2))<1e-14
    assert report["retrievedEgPopulation"]>0.9998
    n=loadmat(DATA/"nominal_trajectories.mat",simplify_cells=True)["nominal"]
    for name,r in n.items():
        for field in ["F","FB","C","Patom","Pphoton","Poutside"]:
            assert np.isfinite(r[field]).all(),(name,field)
            assert r[field].min()>-1e-10 and r[field].max()<1+1e-8,(name,field)
        assert np.max(np.abs(r["Patom"]+r["Pphoton"]-1))<1e-8
        state=loadmat(DATA/f"nominal_{name}_states.mat",simplify_cells=True)["trajectory"]
        assert np.max(np.abs(np.sum(abs(state["psi"])**2,axis=0)-1))<1e-8
    assert np.max(abs(n["dressed"]["FB"]-1))<1e-10
    assert abs(n["attach"]["prepF"]-1)<1e-10
    assert abs(n["attach"]["prepFB"]-report["analyticZ"])<1e-10
    assert n["atomic"]["F"][n["atomic"]["time"]>=20].min()>0.9958
    assert np.ptp(n["atomic"]["FB"][n["atomic"]["time"]>=20])<1e-10
    assert 1-n["atomic"]["prepFB"]<report["residualBound"]
    nominal=n["atomic"]["prepF"]

    old=np.genfromtxt(ROOT/"data/paper/protocol_parameter_robustness_paper_source_data.csv",
                     delimiter=",",names=True,dtype=None,encoding="utf8")
    loading=old[old["map"]=="loading"]
    symmetry_error=0
    for area in np.unique(loading["x"]):
        row=loading[loading["x"]==area]
        symmetry_error=max(symmetry_error,float(np.ptp(row["final_Bell_fidelity"])))
    assert symmetry_error<1e-7
    area_error=np.max(abs(loading["final_Bell_fidelity"]-nominal*np.cos(np.pi*loading["x"]/2)**2))
    assert area_error<1e-7

    b=csv("beta_error_maps.csv");bc=csv("beta_map_grid_checks.csv")
    assert len(b)==21*21
    assert b["norm_error"].max()<1e-8
    beta_grid_error=0
    for a,c in [(0,0),(-.1,-.1),(-.1,.1),(.1,-.1),(.1,.1)]:
        sel=np.isclose(bc["epsilon_beta1"],a)&np.isclose(bc["epsilon_beta2"],c)
        for field in ["F_preparation","F_after_hold"]:
            beta_grid_error=max(beta_grid_error,float(np.ptp(bc[field][sel])))
    assert beta_grid_error<1e-10
    nominal_beta=(np.isclose(b["epsilon_beta1"],0)&np.isclose(b["epsilon_beta2"],0))
    assert abs(b["F_preparation"][nominal_beta][0]-nominal)<1e-9

    matched=csv("matched_k_benchmark.csv")
    assert len(matched)==21*3
    assert matched["norm_error"].max()<1e-8
    assert np.all(matched["minimum_F_hold"]<=matched["F_after_hold"]+1e-12)
    for im,name in [(1,"atomic"),(2,"fixed"),(3,"attach")]:
        row=matched[(matched["method"]==im)&np.isclose(matched["deltaK_fraction"],0)]
        assert abs(row["F_preparation"][0]-n[name]["prepF"])<1e-10
        assert abs(row["F_after_hold"][0]-n[name]["F"][-1])<1e-10
    sweep=csv("carrier_convergence.csv")
    assert np.all(np.diff(sweep["max_Bell_error"])<0)
    assert sweep["norm_error"].max()<1e-8
    result={"status":"passed","nominal_grid_tolerance":1e-10,
            "maximum_beta_grid_difference":beta_grid_error,
            "loading_phase_invariance_error":symmetry_error,
            "loading_area_scaling_error":float(area_error),
            "dressed_tracking_error":float(np.max(abs(n["dressed"]["FB"]-1)))}
    (DATA/"artifact_checks.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))


if __name__=="__main__":main()
