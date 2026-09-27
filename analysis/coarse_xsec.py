#!/usr/bin/env python3
"""Diagnostics for the p(pi+) cross section: per-bin MC statistics, and a
coarse-binned cross section where the acceptance is actually determined.
"""
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(_HERE, "..", "python"))
from calc_xsections import (  # noqa: E402
    read_csv, cut_for_MM, prep_for_ana, hist_data, virtual_photon_flux,
    luminosity, statistical, get_error_bars, theta_bins, w_bins_e99, q2_bins_e99,
)

mc_path, data_path = sys.argv[1], sys.argv[2]
out_csv = sys.argv[3] if len(sys.argv) > 3 else "coarse_xsec.csv"

print("loading...", flush=True)
mc_rec, thrown = read_csv(mc_path)
rec = read_csv(data_path, True)
print(f"  data {len(rec):,}  mc_rec {len(mc_rec):,}  thrown {len(thrown):,}", flush=True)

rec, mc_rec, _ = cut_for_MM(rec, mc_rec, rec.iloc[0:0])
print(f"after MM cut: data {len(rec):,}  mc_rec {len(mc_rec):,}", flush=True)

R = prep_for_ana(rec.copy(), w_bins_e99, q2_bins_e99, theta_bins)
M = prep_for_ana(mc_rec.copy(), w_bins_e99, q2_bins_e99, theta_bins)
T = prep_for_ana(thrown.copy(), w_bins_e99, q2_bins_e99, theta_bins)
print(f"after binning: data {len(R):,}  mc_rec {len(M):,}  thrown {len(T):,}", flush=True)

print("\n=== MC statistics per (W,Q2,cos) bin, summed over phi ===")
print(f"{'W':>13s} {'Q2':>13s} {'cos':>11s} {'mc_rec':>9s} {'thrown':>11s} {'acc':>8s} {'1/sqrt(Nrec)':>12s}")
rows = []
for w in sorted(pd.unique(R.w_bin).dropna()):
    for q2 in sorted(pd.unique(R.q2_bin).dropna()):
        for ct in sorted(pd.unique(R.theta_bin).dropna()):
            m = int(((M.w_bin == w) & (M.q2_bin == q2) & (M.theta_bin == ct)).sum())
            t = int(((T.w_bin == w) & (T.q2_bin == q2) & (T.theta_bin == ct)).sum())
            if t == 0:
                continue
            acc = m / t
            rel = 1 / np.sqrt(m) if m > 0 else np.inf
            rows.append((w, q2, ct, m, t, acc, rel))
            if m < 200:
                print(f"{str(w):>13s} {str(q2):>13s} {str(ct):>11s} {m:9,d} {t:11,d} {acc:8.4f} {rel:12.1%}")
accs = np.array([r[5] for r in rows])
print(f"\nbins total {len(rows)}, median acceptance {np.median(accs):.4f}, "
      f"bins with >=100 mc_rec: {sum(1 for r in rows if r[3] >= 100)}")

# ---- coarse cross section: one W, one Q2 bin, all cos and phi ----
print("\n=== coarse cross sections (all cos, all phi, 1 phi bin) ===")
print(f"{'W bin':>13s} {'Q2 bin':>13s} {'N_data':>9s} {'N_mc':>8s} {'N_thr':>10s} {'acc':>8s} "
      f"{'sigma (ub)':>12s} {'stat%':>7s} {'acc%':>7s} {'tot%':>7s}")
Q_FULL = 4348.46636e-6
lum = luminosity()
out = []
for w in sorted(pd.unique(R.w_bin).dropna()):
    for q2 in sorted(pd.unique(R.q2_bin).dropna()):
        d = R[(R.w_bin == w) & (R.q2_bin == q2) & (R.cut_fid)]
        m = M[(M.w_bin == w) & (M.q2_bin == q2)]
        t = T[(T.w_bin == w) & (T.q2_bin == q2)]
        nd, nm, nt = len(d), len(m), len(t)
        if nd < 20 or nt == 0:
            continue
        acc = nm / nt
        dw, dq = w.right - w.left, q2.right - q2.left
        dcos, dphi = 2.0, 2 * np.pi
        kb = dw * dq * dcos * dphi
        flux = virtual_photon_flux((w.left + w.right) / 2, (q2.left + q2.right) / 2) * lum
        y = (nd / Q_FULL) / (kb * flux * acc)
        st = statistical(np.array([nd]), np.array([0.0]), np.array([kb * flux * acc]))[0]
        acc_err = 1 / np.sqrt(nm)
        tot = np.sqrt(st**2 + (y * acc_err)**2)
        out.append((w, q2, y, st, acc_err))
        print(f"{str(w):>13s} {str(q2):>13s} {nd:9,d} {nm:8,d} {nt:10,d} {acc:8.4f} "
              f"{y:12.4f} {100*st/y:6.1f}% {100*acc_err:6.1f}% {100*tot/y:6.1f}%")

if out:
    ys = np.array([o[2] for o in out])
    print(f"\ncoarse sigma range: {ys.min():.4f} - {ys.max():.4f} ub "
          f"(median {np.median(ys):.4f})")
    pd.DataFrame([{"w_left": o[0].left, "w_right": o[0].right,
                   "q2_left": o[1].left, "q2_right": o[1].right,
                   "sigma_ub": o[2], "stat": o[3], "acc_frac_err": o[4]}
                  for o in out]).to_csv(out_csv, index=False)
    print("wrote", out_csv)
