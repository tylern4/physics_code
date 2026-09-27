#!/usr/bin/env python3
"""Plot the coarse p(pi+) cross section: original vs fixed, and the ratio."""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PRIS = sys.argv[1] if len(sys.argv) > 1 else "coarse_pris.csv"
FIXED = sys.argv[2] if len(sys.argv) > 2 else "coarse_fixed.csv"
OUTDIR = sys.argv[3] if len(sys.argv) > 3 else "xsec_results"

pris = pd.read_csv(PRIS)
fixed = pd.read_csv(FIXED)
m = pris.merge(fixed, on=["w_left", "w_right", "q2_left", "q2_right"], suffixes=("_p", "_f"))
m["ratio"] = m.sigma_ub_f / m.sigma_ub_p
m["tot_err"] = 100 * np.sqrt((m.stat_f / m.sigma_ub_f) ** 2 + m.acc_frac_err_f ** 2)

os.makedirs(OUTDIR, exist_ok=True)
wcent = (m.w_left + m.w_right) / 2
qcent = (m.q2_left + m.q2_right) / 2


def grid(col):
    G = np.full((len(np.unique(qcent)), len(np.unique(wcent))), np.nan)
    Q = {v: i for i, v in enumerate(np.unique(qcent))}
    W = {v: i for i, v in enumerate(np.unique(wcent))}
    for _, r in m.iterrows():
        G[Q[(r.q2_left + r.q2_right) / 2], W[(r.w_left + r.w_right) / 2]] = r[col]
    return G, list(W.keys()), list(Q.keys())


fig, axes = plt.subplots(1, 3, figsize=(19, 5.4))

for ax, col, title in [
    (axes[0], "sigma_ub_p", "original  (1252611)"),
    (axes[1], "sigma_ub_f", "FIXED  (287c056)"),
]:
    G, W, Q = grid(col)
    im = ax.pcolormesh(W, Q, G, shading="nearest", cmap="viridis")
    ax.set_xlabel("W (GeV)")
    ax.set_ylabel("Q$^2$ (GeV$^2$)")
    ax.set_title(f"d$\\sigma$/dW dQ$^2$ ($\\mu$b), {title}", fontsize=10)
    fig.colorbar(im, ax=ax, fraction=0.046, label="$\\mu$b / (GeV$^2$)")

G, W, Q = grid("ratio")
im = axes[2].pcolormesh(W, Q, G, shading="nearest", cmap="coolwarm", vmin=1.0, vmax=1.6)
axes[2].set_xlabel("W (GeV)")
axes[2].set_ylabel("Q$^2$ (GeV$^2$)")
axes[2].set_title("FIXED / original", fontsize=10)
fig.colorbar(im, ax=axes[2], fraction=0.046)
fig.suptitle("p($\\pi^+$) exclusive cross section, E1D, all cos$\\theta_*$ and all $\\phi_*$\n"
             "integrated over cos and phi; acceptance from 300 MC files (18.6M thrown, 325k rec)",
             fontsize=12)
fig.tight_layout(rect=(0, 0, 1, 0.88))
fig.savefig(os.path.join(OUTDIR, "xsec_coarse_map.png"), dpi=100)
plt.close(fig)
print(f"wrote {OUTDIR}/xsec_coarse_map.png")

# W profile at fixed Q2, both results overlaid with error bars
fig, axes = plt.subplots(1, 3, figsize=(17, 4.6), sharey=True)
for ax, (qlo, qhi) in zip(axes, [(1.33, 1.56), (1.56, 1.87), (1.87, 2.23)]):
    s = m[(m.q2_left == qlo) & (m.q2_right == qhi)].sort_values("w_left")
    ax.errorbar(s.w_left, s.sigma_ub_p, yerr=s.stat_p, fmt="-o", ms=3, lw=1.3,
                color="tab:red", label="original")
    ax.errorbar(s.w_left, s.sigma_ub_f, yerr=s.stat_f, fmt="-s", ms=3, lw=1.3, ls="--",
                color="tab:blue", label="FIXED")
    ax.set_xlabel("W (GeV)")
    ax.set_title(f"Q$^2$ = {qlo:.2f}–{qhi:.2f} GeV$^2$", fontsize=11)
    ax.grid(alpha=0.25)
    ax.legend(fontsize=9)
axes[0].set_ylabel("d$\\sigma$/dW dQ$^2$ ($\\mu$b)")
fig.suptitle("p($\\pi^+$) E1D cross section vs W, original vs FIXED "
             "(error bars: data statistics only; acceptance error is 1.5–4%)", fontsize=12)
fig.tight_layout(rect=(0, 0, 1, 0.90))
fig.savefig(os.path.join(OUTDIR, "xsec_coarse_W.png"), dpi=100)
plt.close(fig)
print(f"wrote {OUTDIR}/xsec_coarse_W.png")

m.to_csv(os.path.join(OUTDIR, "coarse_xsec_orig_vs_fixed.csv"), index=False)
print(f"wrote {OUTDIR}/coarse_xsec_orig_vs_fixed.csv ({len(m)} bins)")
print(f"ratio: median {m.ratio.median():.3f}  range {m.ratio.min():.3f}-{m.ratio.max():.3f}")
