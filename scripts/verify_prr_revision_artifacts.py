"""Numerical consistency checks on the saved PRR revision artifacts."""
from pathlib import Path
import json
import numpy as np
from scipy.io import loadmat

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/"data"/"prr_revision"


def csv(name):
    return np.genfromtxt(DATA/name,delimiter=",",names=True,dtype=None,encoding="utf8")


def verify_geometric_appendix():
    prefix="geometric_appendix/"
    a=csv(prefix+"amplitude_rule.csv")
    lam=a["lambda"]
    assert a["eigen_residual"].max()<1e-12
    assert a["reconstruction_error"].max()<1e-12
    assert np.max(abs(a["C_conditional"]-2*lam/(1+lam**2)))<1e-12
    assert np.max(abs(a["C_unconditional"]-a["Z"]*a["C_conditional"]))<1e-12
    r=csv(prefix+"radiation_zero.csv")
    selected=r["delta_k"]<=1e-3
    slopes={}
    for root in ["plus","minus"]:
        for channel,power in [("selected",4),("other",2)]:
            field=f"{channel}_{root}_root"
            slope=float(np.polyfit(np.log(r["delta_k"][selected]),np.log(r[field][selected]),1)[0])
            assert abs(slope-power)<0.01
            slopes[field]=slope
    d=csv(prefix+"dressing_error.csv")
    for n2,coefficient in [(6,0.5),(10,15/17)]:
        row=d[d["n2"]==n2]
        chi=coefficient*(0.9*row["g_over_xi"])**2
        assert np.max(abs(row["Z"]-1/(1+chi)))<1e-12
        assert np.max(abs(row["vacuum_state_infidelity"]-chi/(1+chi)))<1e-12
    m=csv(prefix+"frozen_memory.csv")
    assert m["exact_norm_error"].max()<1e-9
    for T in [20,80,160]:
        row=m[m["T"]==T]
        low=row[np.isclose(row["g_over_xi"],0.01)]["max_atomic_error"][0]
        high=row[np.isclose(row["g_over_xi"],0.02)]["max_atomic_error"][0]
        assert abs(high/low-4)<0.02
    c=csv(prefix+"convergence.csv")
    largest=0
    for T in [20,160]:
        for g in [0.1,0.4]:
            row=c[(c["T"]==T)&np.isclose(c["g_over_xi"],g)]
            assert len(row)==4
            largest=max(largest,float(np.ptp(row["max_atomic_error"])))
    assert largest<1e-9
    return {"geometric_eigen_residual":float(a["eigen_residual"].max()),
            "radiation_zero_powers":slopes,"frozen_memory_convergence_difference":largest}


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
    result.update(verify_geometric_appendix())
    (DATA/"artifact_checks.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))


if __name__=="__main__":main()
