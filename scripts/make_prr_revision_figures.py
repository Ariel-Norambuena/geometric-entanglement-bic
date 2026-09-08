"""Regenerate PRR vector figures from archived and revision source data.

No curve smoothing or fitted values are used. Heatmaps show the computed
grid cells; contours are linear interpolation between those same samples.
Run from any directory after the MATLAB revision scripts finish.
"""
from pathlib import Path
import json
import numpy as np
from scipy.io import loadmat
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import Normalize

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "prr_revision"
OLD = ROOT / "data" / "paper"
OUT = ROOT / "figures" / "prr_revision"
OUT.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({
    "font.family": "STIXGeneral", "mathtext.fontset": "stix",
    "font.size": 10.5, "axes.labelsize": 11, "axes.titlesize": 11,
    "xtick.labelsize": 9.5, "ytick.labelsize": 9.5, "legend.fontsize": 9,
    "axes.linewidth": 0.65, "lines.linewidth": 1.65,
    "xtick.direction": "in", "ytick.direction": "in",
    "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.dpi": 400,
    "axes.spines.top": True, "axes.spines.right": True,
})
BLUE, TEAL, RED, GREY, BLACK = "#0072B2", "#009E73", "#B24772", "#777777", "#222222"
COLORS = {"atomic": BLUE, "dressed": BLACK, "fixed": RED, "attach": TEAL}
NAMES = {"atomic": "Floquet + CD", "dressed": "Dressed ideal",
         "fixed": "Fixed guide", "attach": "Attach after Bell"}
STYLES = {"atomic": "-", "dressed": ":", "fixed": "--", "attach": "-."}


def read_csv(path):
    return np.genfromtxt(path, delimiter=",", names=True, dtype=None, encoding="utf8")


def canvas(rows=2, cols=2, height=5.25):
    fig, axes = plt.subplots(rows, cols, figsize=(7.1, height), layout="constrained")
    fig.get_layout_engine().set(w_pad=0.055, h_pad=0.10, hspace=0.09, wspace=0.07)
    return fig, np.atleast_1d(axes).ravel()


def format_ax(ax, letter, title, xlabel=r"$\xi t$", ylabel=None):
    ax.set_title(f"({letter})  {title}", loc="left", pad=35, weight="normal")
    ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    ax.grid(axis="y", color="#e5e5e5", linewidth=0.45)
    ax.set_axisbelow(True)
    ax.tick_params(length=3, width=0.6)
    ax.margins(x=0)


def legend(ax, ncol=2):
    handles, labels = ax.get_legend_handles_labels()
    if ncol > 2 and sum(map(len, labels)) > 35:
        ncol = 2
    if ncol == 1 and len(labels) > 2:
        ncol = 2
    rows = int(np.ceil(len(labels) / ncol))
    leg = ax.legend(loc="lower left", bbox_to_anchor=(0, 1.015),
                     ncol=ncol, frameon=False, borderaxespad=0,
                     handlelength=2.0, columnspacing=0.85, handletextpad=0.4)
    leg.set_in_layout(False)
    ax.set_title(ax.get_title(loc="left"), loc="left", pad=15+13*rows)
    return leg


def save(fig, name):
    fig.canvas.draw()
    renderer=fig.canvas.get_renderer()
    frame=fig.bbox
    for ax in fig.axes:
        if not hasattr(ax, "_colorbar"):
            assert all(ax.spines[side].get_visible() for side in ["left", "right", "top", "bottom"])
        bounds=[ax.get_tightbbox(renderer)]
        if ax.get_legend() is not None:
            bounds.append(ax.get_legend().get_window_extent(renderer))
        for box in bounds:
            assert box.x0>=frame.x0-1 and box.y0>=frame.y0-1, (name,box)
            assert box.x1<=frame.x1+1 and box.y1<=frame.y1+1, (name,box)
    fig.savefig(OUT / f"{name}.pdf", metadata={"Title": name, "Creator": "Reproducible PRR figure script"})
    fig.savefig(OUT / f"{name}.png")
    plt.close(fig)
    print(name, flush=True)


def two_stage():
    d = read_csv(OLD / "two_stage_bell_bic_protocol_paper_source_data.csv")
    t = d["xi_t"]
    fig, ax = canvas(3, 2, 7.0)
    for field, color, label in [("P_gg_with_CD", GREY, r"$P_{gg}$"),
            ("P_eg_with_CD", BLUE, r"$P_{eg}$"), ("P_ge_with_CD", TEAL, r"$P_{ge}$"),
            ("P_photon_with_CD", RED, r"$P_\gamma$")]:
        ax[0].plot(t, d[field], color=color, label=label)
    format_ax(ax[0], "a", "Populations", ylabel="Population")
    legend(ax[0], 4)
    ax[1].plot(t, d["concurrence_with_CD"], color=BLUE, label="Atomic CD")
    ax[1].plot(t, d["concurrence_without_CD"], color=RED, label="No CD")
    ax[1].plot(t, d["concurrence_target"], ":", color=BLACK, label="Target")
    format_ax(ax[1], "b", "Entanglement", ylabel=r"$\mathcal{C}(t)$")
    legend(ax[1], 3)
    ax[2].plot(t, d["bell_fidelity_with_CD"], color=BLUE, label="Total")
    ax[2].plot(t, d["conditional_bell_fidelity_with_CD"], "--", color=TEAL, label="Conditional")
    ax[2].plot(t, d["bell_fidelity_without_CD"], color=RED, label="No CD")
    format_ax(ax[2], "c", "Bell-state fidelity", ylabel=r"$F_{\Psi^+}(t)$")
    legend(ax[2], 3)
    ax[3].plot(t, d["u1"], color=BLUE, label=r"$u_1$")
    ax[3].plot(t, d["u2"], "--", color=TEAL, label=r"$u_2$")
    format_ax(ax[3], "d", "Floquet coupling path", ylabel=r"$u_j(t)$")
    legend(ax[3])
    ax[4].plot(t, d["Omega_pi_over_xi"], color=BLUE, label="Local resonant pulse")
    format_ax(ax[4], "e", "Loading pulse", ylabel=r"$\Omega_\pi/\xi$")
    ax[4].set_xlim(0, 2)
    legend(ax[4], 1)
    ax[5].plot(t, d["Omega_CD_over_xi"], color=TEAL, label="Atomic exchange")
    format_ax(ax[5], "f", "Counterdiabatic control", ylabel=r"$\dot\theta/\xi$")
    legend(ax[5], 1)
    for a in ax[:4]:
        a.set_ylim(-0.025, 1.03)
    for a in [*ax[:4], ax[5]]:
        a.axvline(2, color=GREY, linestyle=(0, (3, 3)), linewidth=0.9)
        a.set_xlim(0, 22)
    save(fig, "Figure_TwoStage_Bell_BIC")


def validation():
    d = read_csv(OLD / "floquet_effective_validation_paper_source_data.csv")
    f = read_csv(DATA / "carrier_convergence.csv")
    fig, ax = canvas()
    t = d["xi_t"]
    ax[0].plot(t, d["C_effective"], color=BLUE, label="Effective")
    ax[0].plot(t, d["C_lab_slow_frame"], "--", color=TEAL, label="Microscopic")
    format_ax(ax[0], "a", "Concurrence", ylabel=r"$\mathcal{C}(t)$")
    legend(ax[0])
    ax[1].plot(t, d["F_Bell_lab_frame"], color="#aaaaaa", lw=0.65, label="Lab frame")
    ax[1].plot(t, d["F_Bell_effective"], color=BLUE, label="Effective")
    ax[1].plot(t, d["F_Bell_lab_slow_frame"], "--", color=TEAL, label="Slow frame")
    format_ax(ax[1], "b", "Bell fidelity", ylabel=r"$F_{\Psi^+}(t)$")
    legend(ax[1], 3)
    ax[2].plot(t, 1e3*d["abs_delta_bell_fidelity"], color=BLUE, label=r"$10^3|\Delta F_{\Psi^+}|$")
    ax[2].plot(t, 1e3*d["state_infidelity_slow_frame"], "--", color=TEAL, label=r"$10^3(1-F_{\rm state})$")
    format_ax(ax[2], "c", "Resolved micromotion error", ylabel=r"Error $\times\,10^3$")
    legend(ax[2])
    ax[3].semilogy(f["nu_over_xi"], f["max_Bell_error"], "o-", color=BLUE, ms=4, label=r"$\max|\Delta F_{\Psi^+}|$")
    ax[3].semilogy(f["nu_over_xi"], f["max_state_infidelity"], "s--", color=TEAL, ms=4, label=r"$\max(1-F_{\rm state})$")
    format_ax(ax[3], "d", "Carrier-frequency convergence", xlabel=r"$\nu/\xi$", ylabel="Maximum sampled error")
    ax[3].set_xticks([4, 6, 8, 12, 16])
    legend(ax[3])
    for a in ax[:2]: a.set_ylim(-0.02,1.03)
    save(fig, "Figure_Floquet_Effective_Validation")


def dressed_hold():
    n = loadmat(DATA / "nominal_trajectories.mat", simplify_cells=True)["nominal"]
    ret = loadmat(DATA / "retrieval.mat", simplify_cells=True)["retrieval"]
    fig, ax = canvas()
    for mode in ["atomic", "dressed"]:
        r=n[mode]; sel=r["time"]<=20
        ax[0].plot(r["time"][sel], r["FB"][sel], STYLES[mode], color=COLORS[mode], label=NAMES[mode])
    format_ax(ax[0], "a", "Loading the dressed BIC", ylabel=r"$F_B(t)$")
    ax[0].set_ylim(0.999,1.00003)
    ax[0].ticklabel_format(axis="y", style="plain", useOffset=False)
    legend(ax[0], 1)
    for mode in ["atomic", "fixed", "attach", "dressed"]:
        r=n[mode]; sel=r["time"]>=20
        ax[1].plot(r["time"][sel]-20, r["F"][sel], STYLES[mode], color=COLORS[mode], label=NAMES[mode])
    format_ax(ax[1], "b", "Retention with exchange off", xlabel=r"$\xi(t-T)$", ylabel=r"$F_{\Psi^+}(t)$")
    ax[1].set_ylim(0.977,1.001)
    ax[1].ticklabel_format(axis="y", style="plain", useOffset=False)
    legend(ax[1], 1)
    r=n["atomic"]
    ax[2].semilogy(r["time"], np.maximum(r["Pphoton"],1e-14), color=BLUE, label="All guide photons")
    ax[2].semilogy(r["time"], np.maximum(r["Poutside"],1e-14), "--", color=RED, label="Outside sites 0--8")
    ax[2].axvline(20,color=GREY,ls=":",lw=0.9)
    format_ax(ax[2], "c", "Dressing and emitted radiation", ylabel="Photonic population")
    ax[2].set_ylim(1e-8,1e-2)
    legend(ax[2], 1)
    ax[3].plot(ret["time"]-120,ret["egPopulation"],color=BLUE,label="Reversed path + atomic CD")
    ax[3].axhline(1,color=GREY,ls=":",lw=0.9)
    format_ax(ax[3], "d", "Retrieval after the hold", xlabel=r"$\xi(t-T-t_{\rm h})$", ylabel=r"$P_{eg}(t)$")
    ax[3].set_ylim(0.47,1.025)
    legend(ax[3], 1)
    save(fig, "Figure_Dressed_BIC_Retention")


def matched():
    d=read_csv(DATA / "matched_k_benchmark.csv")
    dur=read_csv(DATA / "duration_benchmark.csv")
    fig,ax=canvas()
    for im,mode in enumerate(["atomic","fixed","attach"],1):
        sel=d["method"]==im
        for a,field in [(ax[0],"F_preparation"),(ax[1],"F_after_hold")]:
            a.plot(100*d["deltaK_fraction"][sel],d[field][sel],STYLES[mode],color=COLORS[mode],label=NAMES[mode])
    for i,title in enumerate(["End of preparation","After the same hold"]):
        format_ax(ax[i],chr(97+i),title,xlabel=r"$100\,\delta K/K_0$",ylabel=r"$F_{\Psi^+}$")
        legend(ax[i],1)
    for im,mode in enumerate(["atomic","dressed","fixed","attach"],1):
        sel=dur["method"]==im
        for a,field in [(ax[2],"F_BIC_preparation"),(ax[3],"F_after_hold")]:
            a.plot(dur["T"][sel],dur[field][sel],STYLES[mode],color=COLORS[mode],label=NAMES[mode])
    format_ax(ax[2],"c","BIC loading versus duration",xlabel=r"$\xi T$",ylabel=r"$F_B(T)$")
    format_ax(ax[3],"d","Retention versus duration",xlabel=r"$\xi T$",ylabel=r"$F_{\Psi^+}(T+t_{\rm h})$")
    legend(ax[2],1);legend(ax[3],1)
    for a in ax:
        a.ticklabel_format(axis="y",style="plain",useOffset=False)
    save(fig,"Figure_K_Robustness_Comparison")


def heatmap(ax,x,y,z,letter,title,xlabel,ylabel,norm,contours=(0.99,0.995)):
    mesh=ax.pcolormesh(x,y,z,cmap="cividis",norm=norm,shading="nearest",rasterized=True)
    for level,style in zip(contours,["-","--"]):
        if np.min(z)<level<np.max(z):
            ax.contour(x,y,z,levels=[level],colors=BLACK,linewidths=1,linestyles=style)
    ax.plot(0,0,"o",ms=4.8,mfc="black",mec="white",mew=1)
    ax.set_title(f"({letter})  {title}",loc="left",pad=10)
    ax.set_xlabel(xlabel);ax.set_ylabel(ylabel)
    ax.spines[["top","right"]].set_visible(True)
    return mesh


def effective_maps():
    d=read_csv(OLD / "protocol_parameter_robustness_paper_source_data.csv")
    fig,ax=canvas(height=5.9)
    norm=Normalize(0.9,1)
    specs=[("loading","Loading consistency",100,1,r"$100\,\delta A_\pi/A_\pi$",r"$\delta\varphi_\pi/\pi$"),
        ("floquet_calibration","Effective-coupling errors",100,100,r"$100\,\delta u_1/u_1$",r"$100\,\delta u_2/u_2$"),
        ("timing_cd","Duration and exchange strength",100,100,r"$100\,\delta T/T_0$",r"$100\,\delta\alpha_{\rm CD}$"),
        ("detuning","Atomic detunings",1,1,r"$\Delta_{\rm c}/\xi$",r"$\Delta_{\rm r}/\xi$")]
    for i,(name,title,sx,sy,xlabel,ylabel) in enumerate(specs):
        z=d[d["map"]==name];x=np.unique(z["x"]);y=np.unique(z["y"])
        grid=np.full((len(y),len(x)),np.nan)
        for row in z:
            grid[np.searchsorted(y,row["y"]),np.searchsorted(x,row["x"])]=row["final_Bell_fidelity"]
        assert np.isfinite(grid).all()
        mesh=heatmap(ax[i],sx*x,sy*y,grid,chr(97+i),title,xlabel,ylabel,norm)
    cb=fig.colorbar(mesh,ax=list(ax),fraction=0.033,pad=0.025,shrink=0.92)
    cb.set_label(r"Final unconditional Bell fidelity $F_{\Psi^+}$")
    fig.legend(handles=[Line2D([],[],color=GREY,label=r"$F=0.99$"),
        Line2D([],[],color=GREY,ls="--",label=r"$F=0.995$")],
        loc="outside lower center",ncol=2,frameon=False)
    save(fig,"Figure_Protocol_Parameter_Robustness")


def beta_maps():
    path=DATA / "beta_error_maps.csv"
    if not path.exists():return
    d=read_csv(path);x=np.unique(d["epsilon_beta1"]);y=np.unique(d["epsilon_beta2"])
    fig,ax=canvas(1,2,3.45)
    lower=np.floor(min(d["F_preparation"].min(),d["F_after_hold"].min())*1000)/1000
    upper=np.ceil(max(d["F_preparation"].max(),d["F_after_hold"].max())*1000)/1000
    norm=Normalize(lower,upper)
    for i,(field,title) in enumerate([("F_preparation","After preparation"),("F_after_hold","After the CD-free hold")]):
        grid=np.full((len(y),len(x)),np.nan)
        for row in d:
            grid[np.searchsorted(y,row["epsilon_beta2"]),np.searchsorted(x,row["epsilon_beta1"])]=row[field]
        mesh=heatmap(ax[i],100*x,100*y,grid,chr(97+i),title,r"$100\,\epsilon_{\beta_1}$",r"$100\,\epsilon_{\beta_2}$",norm)
    cb=fig.colorbar(mesh,ax=list(ax),fraction=0.035,pad=0.02)
    cb.set_label(r"$F_{\Psi^+}$")
    save(fig,"Figure_Floquet_Depth_Robustness")


def geometric_appendix():
    folder=DATA / "geometric_appendix"
    r=read_csv(folder / "radiation_zero.csv")
    a=read_csv(folder / "amplitude_rule.csv")
    d=read_csv(folder / "dressing_error.csv")
    m=read_csv(folder / "frozen_memory.csv")
    fig,ax=canvas(height=5.8)
    ax[0].loglog(r["delta_k"],r["selected_plus_root"],color=BLUE,label="Selected channel")
    ax[0].loglog(r["delta_k"],r["other_plus_root"],"--",color=RED,label="Opposite channel")
    format_ax(ax[0],"a","Radiation near the resonant root",xlabel=r"$|k-K_0|$",ylabel=r"$N_c|\mathcal{M}(k)|^2/(g u_0)^2$")
    legend(ax[0],1)
    ax[0].text(0.53,0.27,r"$\propto\delta k^4$",color=BLUE,transform=ax[0].transAxes)
    ax[0].text(0.18,0.76,r"$\propto\delta k^2$",color=RED,transform=ax[0].transAxes)
    aa=a[a["n2"]==10];lam=aa["lambda"]
    ax[1].plot(lam,2*lam/(1+lam**2),color=BLACK,label="Atomic rule")
    ax[1].plot(lam[::3],aa["C_conditional"][::3],"o",color=BLUE,ms=3.8,mfc="white",label="Conditional")
    ax[1].plot(lam,aa["C_unconditional"],"--",color=TEAL,label="Unconditional")
    format_ax(ax[1],"b","Geometric amplitude selection",xlabel=r"$\lambda_F=n_1u_1/(n_2u_2)$",ylabel="Concurrence")
    ax[1].set_ylim(-0.025,1.03)
    legend(ax[1],1)
    for n2,col,ls in [(6,BLUE,"-"),(10,TEAL,"--")]:
        dd=d[d["n2"]==n2]
        ax[2].loglog(dd["g_over_xi"],dd["vacuum_state_infidelity"],ls,color=col,label=rf"$(n_1,n_2)=(6,{n2})$")
    format_ax(ax[2],"c","Accuracy of the vacuum-field state",xlabel=r"$g/\xi$",ylabel=r"$1-|\langle a,0|B\rangle|^2=1-Z$")
    legend(ax[2],1)
    for T,col,ls,marker in [(20,BLUE,"-","o"),(80,TEAL,"--","s"),(160,RED,":","^")]:
        mm=m[m["T"]==T]
        ax[3].loglog(mm["g_over_xi"],mm["max_atomic_error"],ls,color=col,marker=marker,ms=4,mfc="white",label=rf"$\xi T={T}$")
    format_ax(ax[3],"d","Frozen-memory approximation",xlabel=r"$g/\xi$",ylabel=r"$\max_t\|c_{\rm fr}-c_{\rm ex}\|_2$")
    legend(ax[3],3)
    save(fig,"Figure_Geometric_BIC_Approximation")


if __name__=="__main__":
    two_stage();validation();dressed_hold();matched();effective_maps();beta_maps()
    geometric_appendix()
