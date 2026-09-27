# p(π⁺) exclusive cross section — CLAS E1D, deuterium target

Final result from the fixed reconstruction code (`287c056`) on all 480
experimental runs and all 987 simulation files.

## Headline

**dσ/dW dQ² = 0.026–0.576 μb per (W,Q²) bin, median 0.161 μb**, over 212 bins
spanning W = 1.10–1.82 GeV and Q² = 1.1–3.5 GeV². Total error 2.0% median,
4.0% at the 90th percentile, 7.9% worst case.

The clear Δ(1232) peak at W ≈ 1.23 GeV, a second structure near W ≈ 1.5 GeV, and
the expected Q² suppression are all present and well resolved.

## The reconstruction fix raised the cross section by ~20%

Measured on a pinned 300-file MC sample so both builds see identical simulation:

| | original | FIXED | change |
|---|---:|---:|---:|
| data rows in CSV | 96,474 | 121,924 | **+26.4%** |
| MC `thrown` rows | 18,623,696 | 18,623,696 | 0 (bit-identical) |
| MC `mc_rec` rows | 311,114 | 324,576 | **+4.3%** |
| **median σ(FIXED)/σ(original)** | | | **1.200** |

`thrown` comes from the truth branches and is bit-identical, confirming the
difference is confined to the reconstruction path. Range of the ratio:
1.01–1.83, largest (1.6–1.8×) at low W / high Q², ~1.05× above W ≈ 1.7 GeV.

Going from 579 to 987 MC files changed the cross section by only 1.001 median
(acceptance is a ratio, so it is stable) but roughly halved the errors.

## Files

| file | what |
|---|---|
| `coarse_xsec_full_mc.csv` | **final result**: 212 (W,Q²) bins, σ, stat, acceptance error |
| `xsec_coarse_W_full_mc.png` | dσ/dW dQ² vs W at three Q² slices |
| `coarse_xsec_orig_vs_fixed.csv` | pinned-MC original vs fixed, with the ratio |
| `xsec_coarse_map.png` | 2D (W,Q²) maps, original and fixed, plus ratio |
| `xsec_coarse_W.png` | W profiles from the pinned-MC comparison |

The binning is 0.02 GeV in W by ~0.23 GeV² in Q², integrated over all cosθ\* and
all φ.

## MAID2007 is now working, for the correct reaction

`maid_interface` is a Cython extension wrapping the bundled MAID2007 Fortran. It
is required for the φ bin-centering correction, not just for plotting.

Getting it running took four fixes, all of which had silently broken it:

1. **The extension was never built.** `CMakeLists.txt:96` is `#IF(0)` with the
   whole `.pyx` block commented out, and `aao_rad_lund`'s own CMake cannot
   configure under CMake ≥ 4.0 (`FindPythonLibs` was removed). Built via
   `external/aao_rad_lund/build_m/setup_maid.py`, which leaves the submodule
   untouched.
2. **The table filename was a MAID2000 leftover.** `maid_lee` asked for
   `maid2000-NPpi.tbl`, which is not shipped, so the read always failed.
3. **The reaction was wrong.** `maid_2007_Npi` hardcoded `channel_opt=3`, which
   reads `maid07-PNpi.tbl` = p→nπ⁰. It computed π⁰ production, not π⁺.
4. **The tables were in the tarball all along.** `maid07-NPpi.tbl` is the
   p→nπ⁺ channel (MAID writes the target nucleon first, so NPpi is n→pπ⁰, which
   by charge symmetry is p→nπ⁺). There is no p→π⁺p table because MAID2007's
   published channel set is (π⁰p), (π⁰n), (π⁺n), (π⁻p) — p→π⁺p is I₃=+3/2 and
   is not an independent channel.

Verified by the Δ++ peak: σ rises to 16.9 at W = 1.23 GeV and falls to ~1 beyond
W = 1.4, which is the expected p→nπ⁺ shape. The table grid (W 1.08–2.00,
Q² 0–5.0) covers the full analysis range.

The magnitude of the correction is **not** small: factors of ~1.0–1.8 across a
10-bin φ binning. `binCetnerCorrection` now raises rather than silently
returning 1.0, since that fallback was hiding an ~80% effect on some φ bins.
Applied to the φ-resolved output it shifts bin values by 0.915 median
(0.69–1.16). It does not affect the coarse numbers here, which integrate over
all φ.

## Two things still missing

1. **No empty-target subtraction.** That run lives at
   `/work/gothelab/clas6/e1d/empty` and is not in `data/exp`, so the
   target-contamination term is zero and these cross sections are
   **overestimates** by roughly the empty-target contribution. Pass
   `--empty <csv>` once available.
2. **The φ-resolved binning is still not a measurement.** At ~2% acceptance the
   per-bin MC error dominates: over all 14,385 bins the median error/value ratio
   is 304, and **zero** bins have error < value. You need roughly 100
   reconstructed MC events per (W,Q²,cosθ\*,φ) bin for a 10% acceptance error,
   i.e. of order 10–30× more simulation than exists. The coarse binning above is
   the first result that is actually determined.

## Reproducing

```bash
export MOM_CORR=$PWD/data/MomCorr NUM_THREADS=8
./build/csv_maker     -e1d data/exp -o xsec_work/data     # 121,924 rows
./build/csv_maker_mc  -e1d data/sim -o xsec_work/mc       # 987 files, 4.6 GB

# one-time: build the MAID Cython extension
cd external/aao_rad_lund/build_m && ../../../../.venv/bin/python setup_maid.py

cd python && ../.venv/bin/python crossSectionsCSV.py \
    --mc ../xsec_work/mc_e1d.csv --data ../xsec_work/data_e1d.csv \
    --folder ../xsec_work/out_maid --phibins 10
```

`calc_xsections` sets `CLAS_PARMS` to `exe/maid_xsection` on import, which is
how the Fortran locates `spp_tbl/` (via `revinm`). It does not change directory,
so relative paths in the caller keep working.

`xsec_work/diagnose.py` produces the coarse binning and the per-bin MC
statistics; `xsec_work/plot_final.py` makes the final plot;
`xsec_work/plot_xsec.py` makes the original-vs-fixed comparison.
