"""
Block Copolymer Phase Diagram Mapping
SCFT + Leibler theory for PS-b-PMMA equilibrium morphologies

Predicts: DIS -> BCC -> HEX -> Gyroid -> LAM as function of (f_PS, χN)
References:
  - Leibler, Macromolecules 1980
  - Matsen & Bates, Macromolecules 1996
  - Fredrickson & Helfand, J.Chem.Phys 1987
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import json, os

def chi_ps_pmma(T_K):
    return 0.04 + 4.9 / T_K

def F_Leibler(f, x):
    def g(f_, x_):
        return 2.0*(np.exp(-f_*x_) - 1.0 + f_*x_)/(x_**2)
    g1 = g(f, x); g2 = g(1-f, x)
    g12 = 0.5*(g(1.0,x) - g1 - g2)
    W = g1*g2 - g12**2
    return W / (g1 + g2 + 2*g12)

def chi_ODT_MF(f):
    x_arr = np.linspace(0.1, 100.0, 2000)
    F_arr = [F_Leibler(f, x) for x in x_arr]
    idx = np.argmin(F_arr)
    return 0.5 / F_arr[idx], x_arr[idx]

def chi_ODT_fluctuation(f, N):
    chiN_mf, _ = chi_ODT_MF(f)
    if abs(f - 0.5) < 0.02:
        return chiN_mf + 41.022 * N**(-1/3)
    return chiN_mf + 30.0 * N**(-1/3) * (1 - 4*(f-0.5)**2)

def domain_spacing_ssl(N, b, f, chi):
    return 1.03*b*N**(2/3)*chi**(1/6)*(f*(1-f))**(2/3)

def domain_spacing_wsl(N, b):
    Rg = b*np.sqrt(N/6.0)
    return 2*np.pi*Rg/1.95

def estimate_L0(N, T_K=500.0, b=0.68, f=0.5):
    chi = chi_ps_pmma(T_K); chiN = chi*N
    d_ssl = domain_spacing_ssl(N, b, f, chi)
    d_wsl = domain_spacing_wsl(N, b)
    alpha = np.clip((chiN-10.5)/(chiN+10.5), 0, 1)
    return {"N":N,"T_K":T_K,"chi":chi,"chiN":chiN,
            "d_WSL_nm":d_wsl,"d_SSL_nm":d_ssl,
            "d_estimated_nm": alpha*d_ssl+(1-alpha)*d_wsl,
            "regime":"SSL" if chiN>20 else "WSL"}

def compute_full_phase_diagram(N=100):
    f_arr = np.linspace(0.05, 0.95, 181)
    chiN_mf, chiN_fh, q_star = [], [], []
    for f in f_arr:
        cm, xs = chi_ODT_MF(f)
        cf = chi_ODT_fluctuation(f, N)
        chiN_mf.append(cm); chiN_fh.append(cf); q_star.append(xs)
    pb = {
        "BCC_HEX":{"f_L":[0.20,0.23,0.27,0.29,0.30],
                   "f_R":[0.80,0.77,0.73,0.71,0.70],
                   "chiN":[15.0,20.0,30.0,40.0,60.0]},
        "HEX_GYR":{"f_L":[0.26,0.29,0.31,0.33,0.34],
                   "f_R":[0.74,0.71,0.69,0.67,0.66],
                   "chiN":[15.0,20.0,30.0,40.0,60.0]},
        "GYR_LAM":{"f_L":[0.34,0.36,0.38,0.39,0.40],
                   "f_R":[0.66,0.64,0.62,0.61,0.60],
                   "chiN":[15.0,20.0,30.0,40.0,60.0]},
    }
    sym_mf, _ = chi_ODT_MF(0.5)
    return {"f_arr":f_arr.tolist(),"chiN_ODT_MF":chiN_mf,"chiN_ODT_FH":chiN_fh,
            "q_star":q_star,"N":N,"phase_boundaries":pb,
            "symmetric_ODT":{
                "f=0.5_chiN_MF": sym_mf,
                "f=0.5_chiN_FH_N100": chi_ODT_fluctuation(0.5,100),
                "f=0.5_chiN_FH_N200": chi_ODT_fluctuation(0.5,200),
            }}

def plot_phase_diagram(pd_data, output_path="figures/phase_diagram.png"):
    fig, ax = plt.subplots(figsize=(10,7))
    f   = np.array(pd_data["f_arr"])
    mf  = np.array(pd_data["chiN_ODT_MF"])
    fh  = np.array(pd_data["chiN_ODT_FH"])

    ax.plot(f, mf, 'k--', lw=1.5, label="ODT (Mean-Field)")
    ax.plot(f, fh, 'k-',  lw=2.0, label="ODT (Fluctuation-corrected)")

    pb = pd_data["phase_boundaries"]
    styles = [("BCC_HEX","#2E86C1","-"),("HEX_GYR","#28B463","--"),("GYR_LAM","#E74C3C",":")]
    for key, col, ls in styles:
        ax.plot(pb[key]["f_L"], pb[key]["chiN"], color=col, lw=2, ls=ls)
        ax.plot(pb[key]["f_R"], pb[key]["chiN"], color=col, lw=2, ls=ls)

    labels = [
        (0.50,30,"Lamellae (LAM)","#C0392B",13),
        (0.37,38,"Gyroid (GYR)","#B7950B",11),
        (0.29,47,"Cylinders (HEX)","#1E8449",11),
        (0.19,57,"Spheres (BCC)","#1A5276",11),
        (0.63,38,"Gyroid (GYR)","#B7950B",11),
        (0.71,47,"Cylinders (HEX)","#1E8449",11),
        (0.81,57,"Spheres (BCC)","#1A5276",11),
        (0.50,12,"Disordered (DIS)","#808080",13),
    ]
    for xp,yp,txt,col,sz in labels:
        ax.text(xp,yp,txt,ha='center',va='center',fontsize=sz,
                color=col,fontweight='bold',alpha=0.85)

    ax.scatter([0.50],[50.0],color="red",s=180,zorder=10,marker="*",
               label="7nm node target (N≈1000)")
    ax.annotate("7nm node", xy=(0.50,50.0), xytext=(0.62,56),
                fontsize=10, color="darkred",
                arrowprops=dict(arrowstyle="->",color="darkred"))

    ax.set_xlabel("Volume Fraction of PS  ($f_{PS}$)", fontsize=14)
    ax.set_ylabel(r"Interaction Parameter  $\chi N$", fontsize=14)
    ax.set_title("PS-b-PMMA Block Copolymer Phase Diagram\n(Leibler / Matsen-Bates SCFT)",fontsize=14)
    ax.set_xlim(0.05,0.95); ax.set_ylim(8,80)
    ax.legend(loc="upper right",fontsize=10); ax.grid(alpha=0.3)
    os.makedirs(os.path.dirname(output_path),exist_ok=True)
    fig.tight_layout(); fig.savefig(output_path,dpi=200,bbox_inches="tight")
    plt.close(fig); print(f"  Saved: {output_path}")

def plot_domain_spacing(output_path="figures/domain_spacing.png"):
    fig, axes = plt.subplots(1,2,figsize=(13,5))
    N_arr = np.arange(50,1500,10)
    T_arr = [400,450,500,550]
    colors = plt.cm.plasma(np.linspace(0.2,0.85,len(T_arr)))

    ax = axes[0]
    for T,col in zip(T_arr,colors):
        L0 = [estimate_L0(N,T)["d_estimated_nm"] for N in N_arr]
        ax.plot(N_arr,L0,color=col,lw=2,label=f"T={T}K")
    ax.axhline(7,  color="red",   ls="--",lw=1.5,alpha=0.8,label="7nm half-pitch")
    ax.axhline(14, color="orange",ls="--",lw=1.5,alpha=0.8,label="14nm half-pitch")
    ax.axhline(25, color="green", ls="--",lw=1.5,alpha=0.8,label="25nm half-pitch")
    ax.set_xlabel("Total Polymerization Degree  N",fontsize=12)
    ax.set_ylabel("Natural Period  $L_0$  (nm)",fontsize=12)
    ax.set_title("Domain Spacing vs Chain Length",fontsize=13)
    ax.legend(fontsize=9); ax.set_xlim(50,1500); ax.set_ylim(0,50); ax.grid(alpha=0.3)

    ax2 = axes[1]
    f_arr = np.linspace(0.2,0.8,200)
    for N,col in zip([100,200,400,800,1200],
                     plt.cm.viridis(np.linspace(0.1,0.9,5))):
        chi = chi_ps_pmma(500.0)
        L0 = [domain_spacing_ssl(N,0.68,f,chi)*10 for f in f_arr]
        ax2.plot(f_arr,L0,color=col,lw=2,label=f"N={N}")
    ax2.set_xlabel("Volume Fraction  $f_{PS}$",fontsize=12)
    ax2.set_ylabel("Domain Spacing  $L_0$  (Å)",fontsize=12)
    ax2.set_title("L0 vs Composition (T=500K, SSL)",fontsize=13)
    ax2.legend(fontsize=9); ax2.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(output_path,dpi=200,bbox_inches="tight")
    plt.close(fig); print(f"  Saved: {output_path}")

if __name__ == "__main__":
    print("=== Phase Diagram: PS-b-PMMA ===")
    pd = compute_full_phase_diagram(N=100)
    sym = pd["symmetric_ODT"]
    print(f"  ODT (f=0.5) MF:  chiN* = {sym['f=0.5_chiN_MF']:.3f}")
    print(f"  ODT (f=0.5) FH N=100: {sym['f=0.5_chiN_FH_N100']:.3f}")
    print(f"  ODT (f=0.5) FH N=200: {sym['f=0.5_chiN_FH_N200']:.3f}")

    print("\nDomain spacing estimates:")
    spacing_table = []
    for N in [100,200,400,800,1200]:
        res = estimate_L0(N, T_K=500.0)
        print(f"  N={N:4d}  chiN={res['chiN']:.2f}  L0≈{res['d_estimated_nm']:.2f}nm [{res['regime']}]")

    targets = [("7nm node",7.0),("10nm node",14.0),("14nm node",25.0)]
    print("\nSemiconductor process mapping:")
    for label, target in targets:
        for N in range(50,2000):
            res = estimate_L0(N, T_K=500.0)
            if res["d_estimated_nm"] >= target:
                spacing_table.append({"process":label,"target_nm":target,
                    "N_required":N,"chiN":res["chiN"],"L0_nm":res["d_estimated_nm"]})
                print(f"  {label}: N≈{N} chiN≈{res['chiN']:.1f} L0≈{res['d_estimated_nm']:.2f}nm")
                break

    os.makedirs("results",exist_ok=True)
    with open("results/phase_diagram_data.json","w") as f:
        json.dump(pd, f, indent=2)
    with open("results/semiconductor_spacing_table.json","w") as f:
        json.dump(spacing_table, f, indent=2)

    os.makedirs("figures",exist_ok=True)
    plot_phase_diagram(pd)
    plot_domain_spacing()
    print("\n✓ Phase diagram complete.")
