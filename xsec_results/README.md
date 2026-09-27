# p(π⁺) exclusive cross section — CLAS E1D, proton target

Reaction: **e p → e n π⁺**, i.e. MAID's (π⁺, n) channel.

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

## MAID2007 is now building, used as a black box

`maid_interface` is a Cython extension wrapping the bundled MAID2007 Fortran,
and it is required for the φ bin-centering correction, not just for plotting.
No Fortran or Cython source was changed to get it working — the only edit is to
the Cython block of `aao_rad_lund/CMakeLists.txt`.

The module was never built anywhere. This repo's `CMakeLists.txt:96` is `#IF(0)`
with the whole `.pyx` block commented out, and the submodule's own CMake could
not configure under CMake 4.3: it called `find_package(PythonLibs)` and
`FindPythonExtensions.cmake`, both removed in CMake 4, which left
`PYTHON_INCLUDE_DIRS` empty and aborted with "list GET given empty list".
Replacing those with `find_package(Python3 COMPONENTS Interpreter
Development.Module NumPy)` fixes it. The Fortran runtime also has to be linked
explicitly, since a C++-linked MODULE library does not pull in `libgfortran` and
the module otherwise builds but fails to `dlopen`.

```bash
cmake -S external/aao_rad_lund -B external/aao_rad_lund/build \
      -DPython3_EXECUTABLE=$PWD/.venv/bin/python \
      -DCYTHON_EXECUTABLE=$PWD/.venv/bin/cython
cmake --build external/aao_rad_lund/build
```

`CYTHON_EXECUTABLE` must be passed explicitly when Cython is in a venv, or
`find_package(Cython)` reports not-found and the `.pyx` files are skipped with no
error at all.

Verified unchanged: `maid_2007_Npi(4.81726, 1.23, 0.5, 0.0, 90.0)` = 18.71912
before and after the build fix.

It reads `spp_tbl/maid07-PNpi.tbl`, which is the correct table for this
reaction. MAID writes the target nucleon first, so `PN` is proton target with a
neutron in the final state, and with `m_pi = m_pip` assigned at
`maid_lee.f90:20` that is the π⁺ channel. There is no `maid07-PPpi.tbl` and
there never was one — MAID2007's published channels are (π⁰p), (π⁰n), (π⁺n),
(π⁻p), and a proton target can only reach (π⁺n) and (π⁰p) by charge and I₃
conservation. The table's grid (W 1.08–2.00, Q² 0–5.0) covers the whole analysis
range, and the Δ++ shape is right: σ peaks at W = 1.23 GeV and falls away above
W = 1.4.

## What the correction actually does

The bin-centering correction is **not** a small effect. Switching it on, over
all 14,385 φ-resolved bins, the ratio to the uncorrected value has median 0.978
but spans **0.15 to 2.43** — it is a large per-bin reshuffling, not a uniform
shift. The spread is roughly flat in Q² (median ratio 0.99 at Q² = 1.1 down to
0.94 at Q² = 2.7).

Two things about the interface are worth knowing before trusting it:

- `maid_lee.f90` reads its parameter table **once per process** (the `n_call`
  guard at line 29), so a second call with a different `channel_opt` silently
  reuses the first table. Channel comparisons must be run in separate processes.
- `maid_dsigma_all` cannot be called from ctypes: the `extern "C"` block in
  `maidinterface.cpp` is never closed (stray brace at line 88), so that symbol is
  C++-mangled while `maid_dsigma`, defined after the brace, is correctly
  C-linkage. Calling the raw Fortran entry point `maid_` with `byref` scalars
  works. Also, `helicity=0` returns the *unpolarised* cross section
  (`xsection.f90:58`, `sigma_0 = sigma_u + e_hel*sigma_p`) and forces
  `asym_p` to 0 (`maid_lee.f90:370`), so the beam-spin asymmetry is not
  available from this interface.

None of this affects the coarse numbers below, which integrate over all φ.

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

# one-time: build the MAID Cython extension (see the section above)
cmake -S external/aao_rad_lund -B external/aao_rad_lund/build \
      -DPython3_EXECUTABLE=$PWD/.venv/bin/python \
      -DCYTHON_EXECUTABLE=$PWD/.venv/bin/cython
cmake --build external/aao_rad_lund/build

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
