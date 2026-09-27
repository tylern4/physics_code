# Handoff: CLAS E1D p(π⁺) cross section — container + full-chain script

Container work is unfinished: ROOT has to be built from source, and the build
was never taken past configuring. Target for the next session is **NERSC
Perlmutter (x86_64 Linux)**. Nothing in the Dockerfile or the chain script is
specific to any one machine, so the whole thing should carry over unchanged.

---

## 1. What the analysis is

Reaction: **e p → e n π⁺** (proton target, neutron in the final state). This is
MAID's **(π⁺, n)** channel.

Full chain:

```
ROOT files (exp + sim)  ──►  csv_maker / csv_maker_mc  ──►  CSVs
                       ──►  crossSectionsCSV.py         ──►  φ-resolved xsec
                       ──►  analysis/coarse_xsec.py     ──►  coarse xsec + plots
```

Physics result already obtained (do not redo, see `xsec_results/`):

| quantity | value |
|---|---|
| cross section, integrated over all cosθ\* and φ | 0.026 – 0.576 μb per (W,Q²) bin |
| median | 0.161 μb |
| number of bins | 212 |
| total error | 2.0 % median, 4.0 % at 90th pct, 7.9 % worst |
| MC sample used | all 987 simulation files (4.65 GB CSV) |
| Δ(1232) peak | W ≈ 1.23 GeV |
| second structure | W ≈ 1.5 GeV |

Shape is correct: Δ(1232) at W≈1.23, a second bump at W≈1.5, and the expected
Q² suppression.

---

## 2. Git state

Branch `perf/optimize-hot-path`. Working tree is **clean except for the
uncommitted container files** listed in §6.

```
c4d0c4c  Correct the reaction in the results writeup: proton target, e p -> e n pi+
fa29802  Bump aao_rad_lund: Cython modules build again under CMake 4
c0ab849  Build the MAID Cython extension through the submodule's own CMake
3b88bd3  Restore the MAID2007 phi bin-centering correction
1037aac  Make the cross-section pipeline runnable again
287c056  Bound the particle index by _gpart in the hit-number accessors   <- the main fix
```

Submodule `external/aao_rad_lund` is at `69f2882` (pushed to
`github.com/tylern4/aao_rad_lund`). Its only change from upstream is the Cython
block of `CMakeLists.txt` plus a `.gitignore` line for `build/`.

---

## 3. The main physics fix (already done, keep it)

Root cause: commits `d5c5885` (DC/EC/SC, 2018-12-12) and `010508b` (CC,
2019-01-15) changed H10 hit-number indexing to hit-number indirection **without
updating the bounds checks**. Fixed in `287c056` by using `i < _gpart` and
`(unsigned)(hit - 1) < (unsigned)hit_part` across 46 accessors in
`src/branches.cpp`.

Impact, measured on a pinned 300-file MC sample so both builds see identical
simulation:

| | original | fixed | change |
|---|---:|---:|---:|
| data rows in CSV | 96,474 | 121,924 | **+26.4 %** |
| MC `thrown` rows | 18,623,696 | 18,623,696 | 0 (bit-identical) |
| MC `mc_rec` rows | 311,114 | 324,576 | **+4.3 %** |
| median σ(fixed)/σ(original) | | | **1.200** (range 1.01–1.83) |

`thrown` comes from the truth branches and is bit-identical, which confirms the
change is confined to the reconstruction path. The correction is largest
(1.6–1.8×) at low W / high Q², ~1.05× above W ≈ 1.7 GeV.

Verified three full-corpus runs produce all 1,524 histograms bit-identically.

Build flags are plain `-O3`, deliberately **no** `-Ofast` and **no**
`-ffp-contract=fast`; the reasoning is recorded in a long comment in
`CMakeLists.txt` (NaN is used as a sentinel and `-Ofast` makes the compiler
constant-fold `isnan` to false; `-ffp-contract=fast` changes the last digits of
theta and mm2).

---

## 4. Things that cost time — do not re-derive these

### 4.1 MAID tables and which channel is the right one

`maid_2007_Npi` is a **black box** — do not modify the Fortran or the Cython.
It hardcodes `channel_opt=3`, which reads `spp_tbl/maid07-PNpi.tbl`.

* MAID writes the target nucleon first, so `PN` = proton target, neutron final
  state. `maid_lee.f90:20` assigns `m_pi = m_pip` to it (only channel 1, `PPpi`,
  gets `m_pi0`). So `PNpi` **is** the π⁺ channel and the original `channel_opt=3`
  is correct.
* A proton target can only reach (π⁺n) and (π⁰p) by charge and I₃ conservation.
  MAID2007's published channels are (π⁰p), (π⁰n), (π⁺n), (π⁻p), so the set is
  complete.
* **There is no `maid07-PPpi.tbl` and there never was one.** An earlier session
  wrongly searched for a "p→π⁺p" table. It does not exist.
* The only MAID2007 tables distributed are 15 files, neutron- and π⁰n-target
  ones. Fetched in the Dockerfile from
  `https://github.com/tylern4/onepigen/raw/7ad5c99aa5b526dd73311500d680458efb793494/spp_tbl.tar.gz`
  (~118 MB, gitignored, needed at run time).
* Verified: σ peaks at W = 1.23 GeV (the Δ++), table grid W 1.08–2.00,
  Q² 0–5.0 covers the analysis range.
  Reference value: `maid_2007_Npi(4.81726, 1.23, 0.5, 0.0, 90.0) == 18.7191162109375`
  (asserted in the Dockerfile as a smoke test).

### 4.2 The φ bin-centering correction is NOT small

Switching it on changes bin values by a median factor 0.978 but spanning
**0.15 to 2.43** over all 14,385 φ-resolved bins. It is a large per-bin
reshuffling, not a uniform shift. The old code silently returned 1.0 when MAID
was missing, which hid this. `binCetnerCorrection` keeps the 1.0 fallback (its
original signature) but now warns with the size of the effect.

The coarse binning is unaffected — it integrates over all φ.

### 4.3 MAID interface traps

* `maid_lee.f90:29` reads its parameter table **once per process**
  (`if (n_call.ne.0) goto 100`). A second call with a different `channel_opt`
  silently reuses the first table. Channel comparisons need separate processes.
* `maid_dsigma_all` **cannot** be called from ctypes: the `extern "C"` block in
  `maidinterface.cpp` is never closed (stray brace at line 88), so that symbol is
  C++-mangled, while `maid_dsigma` (defined after the brace) is correctly
  C-linkage. The raw Fortran entry point `maid_` works with `byref` scalars.
* `helicity=0` returns the **unpolarised** cross section, not the sum over the
  two helicities (`xsection.f90:58`, `sigma_0 = sigma_u + e_hel*sigma_p`), and
  forces `asym_p = 0` (`maid_lee.f90:370`). Beam-spin asymmetry is therefore not
  available from this interface.
* Tables resolve through `revinm()` using `$CLAS_PARMS`. `calc_xsections` sets it
  to the directory containing `spp_tbl/`; it deliberately does **not** `chdir`,
  because that silently breaks relative paths in the caller (this bit me once).

### 4.4 CSV column conventions (easy to get wrong)

Columns are `electron_sector, w, q2, theta, phi, mm2, cut_fid, helicty, type`.
There is **no header line** in the CSVs.

* `theta` is an **angle in radians** (`_temp->Theta()` in
  `src/reaction.cpp:459`), despite the histogram ranges in
  `include/histogram.hpp` being −1..1 and the Python bin labels saying
  `cos(theta)`. `prep_for_ana` does `np.cos(dataframe.theta)` before binning, so
  the code is self-consistent. `cos(theta)` does land in [−1, 1] for real data
  (checked: −0.9623 .. 1.0000).
* `phi` is in **radians**, 0 .. 2π.
* MC CSV has both `mc_rec` and `thrown` rows in one file, distinguished by the
  `type` column.

### 4.5 Environment variables

* `MOM_CORR` — **required**, the analysis calls `exit(99)` without it. Points at
  the directory holding `angles_*.data` etc. (`data/MomCorr`, 30 files).
* `NUM_THREADS` — default 8.
* `CLAS_PARMS` — set by `calc_xsections`, see §4.3.

### 4.6 Other stale-Python fixes now in place (`1037aac`)

1. `from maid_interface import ...` was a hard import of a module that could not
   be built (§4.5 / §5).
2. `cut_for_MM` returned 2 frames instead of 3 when the empty-target input was
   empty, so the caller failed to unpack. Always 3 now.
3. `pd.Index.sort_values(pd.unique(df.w_bin))` raises under pandas 3
   (`Categorical` has no `is_monotonic_increasing`). Replaced with `sorted()`,
   which also preserves the `Interval` objects the loop needs.
4. `--empty` was `required=True`, making the analysis unrunnable without the
   empty-target run.

---

## 5. Why the Cython module was not being built

It was never built anywhere:

* This repo's `CMakeLists.txt:96` is `#IF(0)` with the whole `.pyx` block
  commented out, so the Cython modules are dead code and the physics compiles
  from the `.cpp` files.
* `aao_rad_lund`'s own CMake could not configure under CMake 4.3: it called
  `find_package(PythonLibs)` and `FindPythonExtensions.cmake`, both **removed in
  CMake 4** (policy CMP0148). That left `PYTHON_INCLUDE_DIRS` empty and
  `FindPythonExtensions.cmake` died with "list GET given empty list".

Fix (committed in the submodule, `69f2882`):

* `find_package(Python3 COMPONENTS Interpreter Development.Module NumPy)`
  replaces the four old calls. `Development.Module` is the modern equivalent of
  `python_extension_module()`: it gives the include dirs and links no libpython.
* The **gfortran runtime must be linked explicitly**. `aao_libs` is mostly
  Fortran, but the generated module is a C++-linked MODULE library, so it builds
  fine and then fails to `dlopen` on undefined `__gfortran_*`.
* `CXX_STANDARD 17` on the generated C++.
* `LIBRARY_OUTPUT_DIRECTORY` so the module appears as `<name>.so` in the build
  dir.

**Trap:** if Cython lives in a venv and is not on `PATH`, `find_package(Cython)`
reports not-found and the `.pyx` files are skipped **with no error at all**. Always
pass `-DCYTHON_EXECUTABLE=<venv>/bin/cython`.

Build:

```bash
cmake -S external/aao_rad_lund -B external/aao_rad_lund/build \
      -DCMAKE_BUILD_TYPE=Release \
      -DPython3_EXECUTABLE=$PWD/.venv/bin/python \
      -DCYTHON_EXECUTABLE=$PWD/.venv/bin/cython
cmake --build external/aao_rad_lund/build
```

---

## 6. Container work — current state

**Files written but NOT yet committed:**

| file | state |
|---|---|
| `Dockerfile` | written, build not yet completing |
| `.dockerignore` | written, uncommitted |
| `scripts/run_full_chain.sh` | written, syntax-checked, never run end-to-end |
| `analysis/coarse_xsec.py` | works on host, verified against known output |
| `analysis/plot_coarse_xsec.py` | works on host, verified |
| `analysis/plot_orig_vs_fixed.py` | works on host, verified |

The three `analysis/*.py` scripts were moved out of the gitignored `xsec_work/`
and their hardcoded `xsec_work/...` paths parameterised, so the chain is
reproducible from git. Verified: `coarse_xsec.py` reproduces
`0.0259 - 0.5756 ub (median 0.1612)` on the 987-file MC.

### 6.1 Dockerfile design

Three stages, all `FROM docker.io/library/ubuntu:24.04` so it is architecture
independent (amd64 and arm64):

1. `rootbuild` — builds **ROOT from source** → `/opt/root`.
2. `build` — venv, `csv_maker`/`csv_maker_mc`, `maid_interface.so`.
3. `runtime` — executables on PATH, `maid_interface.so`, `spp_tbl`, `python/`,
   `analysis/`, `run_full_chain.sh`. Runs as uid/gid 1001.

**ROOT is built from source so one Dockerfile covers both x86_64 and arm64.**
That matters because the prebuilt distributions are not: Ubuntu 24.04 packages
ROOT for amd64 only, and the `rootproject/root` images are amd64-only. Only the
components the code needs are enabled:

* needed: Core, Maths, Physics (TLorentzVector), Hist (TH1D/TH2D), RIO, Net,
  Grafpad (TCanvas/TGraph are referenced)
* disabled: **Cling (`cxx`/`for`) and PyROOT (`python`)** — the two expensive
  parts — plus X11/GL/Java/Qt/Eve/GridFTP/SQL/GiRaFu/DCache.

`ROOT_VERSION=6.34.00` is a build ARG. Note ROOT tags releases with **dashes**
(`v6-34-00`), the Dockerfile translates.

### 6.2 Build failures hit so far, and fixes

1. **`v6.34.00` tag 404** — ROOT tags are `v6-34-00`, not `v6.34.00`. Fixed with
   a `tr '.' '-'` translation.
2. **`curl: (92) HTTP/2 stream ... was not closed cleanly`** fetching the ~90 MB
   ROOT tarball. Fixed with `--http1.1 --retry 5 --retry-delay 3
   --retry-connrefused` (applied to both the ROOT and spp_tbl downloads).
3. **`-DTEST=OFF` breaks ROOT's configure** — ROOT's CMakeLists uses the bare
   word `TEST` as a condition, and a cache variable of that name shadows CMake's
   `if(TEST)` command form. Configure dies at `CMakeLists.txt:289` with
   `if given arguments: "(" "var" "MATCHES" ... Unknown arguments specified`.
   **Remove `-DTEST=OFF`.** This was just fixed in the Dockerfile and is the next
   thing to retry.

**Not yet verified:** that ROOT actually compiles with the trimmed flag set, that
the physics code builds against it, and that `maid_interface.so` links. Those
were never reached.

### 6.3 Notes for running the image on NERSC

Nothing in `Dockerfile` or `scripts/run_full_chain.sh` is tied to a particular
machine, so nothing needs changing to move to Perlmutter. What does need a
decision there:

* **Container runtime.** NERSC's native route is Apptainer, but this image is
  built for OCI runtimes (`podman`/`docker`). Either use `podman build` +
  `podman run` on a login/compute node, or convert for Apptainer with
  `apptainer build e1d.simg docker-daemon://e1d:latest`. The chain script is
  runtime-agnostic.
* **ROOT build time.** Compiling ROOT with Cling and PyROOT disabled is much
  cheaper than a full build, but it is still the long pole. On a login node
  build the image once and push it to a registry, or build in a Slurm job with
  many cores. The ROOT stage needs a few GB of scratch.
* **`--user 1001:1001`.** The image defaults to uid/gid 1001 so files written to
  a shared filesystem keep the right ownership. On NERSC, use your own uid, or
  run as your NERSC user, otherwise output on a shared or project filesystem
  will be unowned. Override with `podman run -u $(id -u):$(id -g)` or
  `--build-arg RUN_UID=... --build-arg RUN_GID=...`.
* **Where the data lives.** See §7. The simulation data must be reachable from
  wherever the container runs; on Perlmutter that is a question of whether the
  share is mounted on the compute nodes and how (CIFS/NFS home, or a
  project-side mount). Confirm the path before wiring it up — a bind mount of a
  path that does not exist can silently produce an empty directory rather than
  an error, which is easy to miss.

---

## 7. Simulation data share

The full simulation sample is on an NFS export, currently observed at
`192.168.86.49:/DangerZone/Danger`, with the files under
`scratch/npip/npip/`. It has been mounted under a few different paths; the
container-visible one is what matters, and it must be established **before**
launching the container (the container cannot mount NFS itself without extra
privileges). Treat it as a setup step and say so out loud rather than doing it
silently.

Contents: `scratch/npip/npip/` holds **55,369 ROOT files of ~26 MB each
(~1.4 TB)**, across many runs (`e1d_osg_*`, `e1d_osg_highq2_*`, ~2000–5000 files
each).

### 7.1 The scale problem — read this before running anything big

Measured rate: the MC CSV is **~4.7 MB per simulation file** (987 files →
4.65 GB). So the full sample needs **~260 GB of intermediate CSV**.

That will not fit in a default workspace, so check free space before starting.
Reading 1.4 TB over NFS is also the dominant wall-clock cost.

The user decided at the end of the last session: **verify the chain on a small
cut; do not run the big job.** Run the full sample on the NERSC machine, or
chunk it and keep only histograms (that needs analysis changes, so it stops
being a black box).

Statistical context for deciding how much MC is worth having: 987 files already
give 2.0 % median error on the coarse result, and the errors scale as
1/√N. The reason to want more MC is the φ-resolved binning, which needs of order
10–30× more than 987 files to become usable. A few thousand files is a
reasonable compromise; 55 k buys very little for the coarse result.

---

## 8. Next steps, in order

1. **Finish the image build.** Retry `podman build -t e1d:latest .` now that
   `-DTEST=OFF` is removed. Expect the ROOT compile to be the long pole (30–90 min
   depending on cores). If the trimmed flag set fails to configure, bisect by
   dropping the `-Dxxx=OFF` group; the ones most likely to be wrong are the
   storage plugins (`gfal`, `dcache`) and `-Drpath=OFF`. The build needs several
   GB of scratch.
2. **Verify the whole chain inside the container on a small cut**, e.g. 200
   files, against the mounted simulation and the local experimental data.
   `$SIMMOUNT` is wherever the share ended up being mounted on the host running
   the container (§7):
   ```bash
   SIMMOUNT=/path/to/mounted/share     # contains scratch/npip/npip
   podman run --rm \
     -v "$PWD/out:/work/out" \
     -v "$PWD/data:/data:ro" \
     -v "$SIMMOUNT:/scratch" \
     -u "$(id -u):$(id -g)" e1d:latest \
     /opt/e1d/bin/run_full_chain.sh \
       --exp /data/exp --sim /scratch/scratch/npip/npip \
       --momcorr /data/MomCorr --out /work/out --cut 200
   ```
   Sanity check: the coarse output should be in the right ballpark
   (≈0.16 μb median) with much larger errors than the 987-file result, and the
   image's own build-time assertion
   `maid_2007_Npi(4.81726,1.23,0.5,0.0,90.0) == 18.7191162109375` must pass.
3. **Then run the full sample** with `--cut` removed, once there is disk for the
   ~260 GB intermediate CSV (§7.1).
4. **Commit** `Dockerfile`, `.dockerignore`, `scripts/run_full_chain.sh`,
   `analysis/`, and update `xsec_results/README.md` with the container usage.
5. Optional / still open: the empty-target run is still missing, so all cross
   sections remain **overestimates** by the target-contamination contribution.
   Pass `--empty <csv>` when it is available.

---

## 9. Reference commands outside the container

Not needed to run the chain, but useful as a baseline for the expected numbers
if you want to compare a container run against a direct run.

```bash
export MOM_CORR=$PWD/data/MomCorr NUM_THREADS=8
./build/csv_maker    -e1d data/exp -o xsec_work/data     # 121,924 rows
./build/csv_maker_mc -e1d data/sim -o xsec_work/mc       # 987 files -> 4.65 GB

.venv/bin/python python/crossSectionsCSV.py \
    --mc xsec_work/mc_e1d.csv --data xsec_work/data_e1d.csv \
    --folder xsec_work/out --phibins 10

.venv/bin/python analysis/coarse_xsec.py \
    xsec_work/mc_e1d.csv xsec_work/data_e1d.csv xsec_work/coarse_xsec.csv
.venv/bin/python analysis/plot_coarse_xsec.py xsec_work/coarse_xsec.csv xsec_results
```

`./build/` must be populated first (top-level `cmake` configure + build against a
ROOT, see the flags in the Dockerfile's `build` stage). The venv is Python
**3.12** — pin it; a newer default was not tested.

Required Python packages: `numpy pandas pyarrow boost-histogram scipy lmfit
matplotlib tqdm cython setuptools`. The chain does **not** need `uproot`, `ROOT`
(PyROOT), `h10` or `nicks_plot_utils` — those appear only in auxiliary scripts
(`python/acceptance.py`, `find_lowest_q2.py`, `vertex_corrections.py`).

---

## 10. Things not to get wrong

* `data/` is the user's data — do not delete it. `xsec_work/` is my own scratch
  and is regenerable; I deleted a 4.3 GB `mc_e1d.csv` from it to free disk.
* ROOT files are read via `glob(pattern + "/*")`, so pass a **directory**
  containing the files, not a glob. It globs all entries, not just `*.root`, so
  the directory must not contain subdirectories.
* `find_package(ROOT REQUIRED COMPONENTS RIO Net)` in the top-level
  `CMakeLists.txt` reads `$ENV{ROOTSYS}` — the Dockerfile sets `ROOTSYS=/opt/root`
  and puts `/opt/root/bin` on `PATH`.
* The physics `CMakeLists.txt` links `${OpenMP_CXX_LIB_NAMES}`, which is
  generally empty; the code uses `std::thread`/`std::async`, not OpenMP, so
  `-fopenmp` is not required.
