"""
Final coarse p(pi+) cross section from the full simulation sample.

Reads the coarse cross section CSV (produced by coarse_xsec.py, one W x one Q2 bin,
integrated over all cos(theta*) and all phi) and writes the deliverable CSV
plus a W-profile plot.

Usage: plot_coarse_xsec.py [coarse_xsec.csv] [outdir]

MAID is not involved here: the bin-centering correction only matters once the
data are resolved in phi, and this binning integrates over the full phi range.
"""

import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

IN = sys.argv[1] if len(sys.argv) > 1 else "coarse_xsec.csv"
OUTDIR = sys.argv[2] if len(sys.argv) > 2 else "xsec_results"
os.makedirs(OUTDIR, exist_ok=True)
d = pd.read_csv(IN)
d["w"] = (d.w_left + d.w_right) / 2.0
d["q2"] = (d.q2_left + d.q2_right) / 2.0
d["tot_err"] = 100 * np.sqrt(d.stat**2 + d.acc_frac_err**2)

d.to_csv(os.path.join(OUTDIR, "coarse_xsec_full_mc.csv"), index=False)
print(f"wrote {OUTDIR}/coarse_xsec_full_mc.csv ({len(d)} bins)")

# Group the Q2 axis into a few representative slices.
slices = [(1.1, 1.35), (1.6, 1.9), (2.3, 2.7)]
labels = [r"$1.1<Q^2<1.35$", r"$1.6<Q^2<1.9$", r"$2.3<Q^2<2.7$"]
colors = ["#1f77b4", "#ff7f0e", "#2ca02c"]

fig, ax = plt.subplots(figsize=(7.2, 5.0))
for (lo, hi), lab, col in zip(slices, labels, colors):
    s = d[(d.q2 >= lo) & (d.q2 < hi)].sort_values("w")
    if s.empty:
        continue
    ax.errorbar(
        s.w, s.sigma_ub, yerr=s.tot_err / 100 * s.sigma_ub,
        marker="o", ms=3.5, lw=1.2, capsize=2, color=col, label=lab,
    )

ax.set_xlabel(r"$W$ (GeV)")
ax.set_ylabel(r"$d\sigma/(dW\,dQ^2)$ ($\mu$b / GeV$^3$)")
ax.set_title("CLAS E1D  $e\\,p \\to e\\,n\\,\\pi^+$, proton target\n"
             "full MC sample, integrated over all $\\cos\\theta^*$ and $\\phi$")
ax.legend(title=r"$Q^2$ (GeV$^2$)", frameon=False)
ax.grid(alpha=0.25)
fig.tight_layout()
fig.savefig(os.path.join(OUTDIR, "xsec_coarse_W_full_mc.png"), dpi=150)
print(f"wrote {OUTDIR}/xsec_coarse_W_full_mc.png")
