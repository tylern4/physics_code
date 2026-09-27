#!/usr/bin/env bash
#
# Run the full p(pi+) cross-section chain:
#
#   ROOT files ---> csv_maker / csv_maker_mc ---> CSVs
#              ---> crossSectionsCSV.py         ---> phi-resolved cross sections
#              ---> analysis/coarse_xsec.py     ---> coarse cross section + MC stats
#              ---> plots
#
# Designed to run inside the container image built from ./Dockerfile, and also
# works on a host that already has the executables and the venv built.
#
#   podman run --rm -v "$PWD/out:/work/out" -v /path/to/data:/data:ro \
#       -u 1001:1001 e1d:latest \
#       /opt/e1d/bin/run_full_chain.sh --exp /data/exp --sim /data/sim \
#                                     --momcorr /data/MomCorr --out /work/out
#
# Every stage is skippable and restartable: a stage whose output already exists
# is skipped unless --force is given.
#
set -euo pipefail

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
EXP=""
SIM=""
MOMCORR=""
OUT="./out"
THREADS="${NUM_THREADS:-8}"
PHIBINS=10
# Limit the number of input files per stream. 0 means no limit. Useful for a
# quick end-to-end check before committing to a multi-hour run over the full
# simulation sample.
CUT=0
FORCE=0
SKIP_XSEC=0
ONLY_EXPERIMENTAL=0

# The phi-resolved binning is 10 (W,Q2) x 10 cos(theta*) x 10 phi bins and is
# not statistically determined at the available Monte Carlo statistics: the
# median error/value ratio is ~300 and no bin has error < value. It is produced
# because it is part of the chain, but analysis/coarse_xsec.py is the output
# that should be quoted.

usage() {
  cat <<'EOF'
Usage: run_full_chain.sh [options]

Input/output:
  --exp DIR        directory of experimental ROOT files (csv_maker)
  --sim DIR        directory of simulation ROOT files (csv_maker_mc)
  --momcorr DIR    directory of momentum-correction tables (data/MomCorr).
                   Required: the analysis exits 99 without it.
  --out DIR        output directory for CSVs, cross sections and plots
                   (default: ./out)

Control:
  --threads N      NUM_THREADS for csv_maker (default: $NUM_THREADS or 8)
  --cut N          only read the first N files per stream (0 = all).
                   Use a small number to smoke-test the whole chain.
  --phibins N      phi bins for the resolved cross section (default: 10)
  --force          recompute stages whose output already exists
  --skip-xsec      skip crossSectionsCSV.py (coarse binning does not need it)
  --exp-only       only run the experimental stream
  -h, --help       this message

Examples:
  # smoke test: 20 files each, no resolved cross section
  run_full_chain.sh --exp data/exp --sim data/sim --momcorr data/MomCorr \
                    --out /tmp/smoke --cut 20 --skip-xsec

  # full simulation sample over the NFS share
  run_full_chain.sh --exp data/exp --sim /scratch/scratch/npip \
                    --momcorr data/MomCorr --out /work/out
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --exp)       EXP="$2"; shift 2 ;;
    --sim)       SIM="$2"; shift 2 ;;
    --momcorr)   MOMCORR="$2"; shift 2 ;;
    --out)       OUT="$2"; shift 2 ;;
    --threads)   THREADS="$2"; shift 2 ;;
    --cut)       CUT="$2"; shift 2 ;;
    --phibins)   PHIBINS="$2"; shift 2 ;;
    --force)     FORCE=1; shift ;;
    --skip-xsec) SKIP_XSEC=1; shift ;;
    --exp-only)  ONLY_EXPERIMENTAL=1; shift ;;
    -h|--help)   usage; exit 0 ;;
    *) echo "unknown option: $1" >&2; usage; exit 2 ;;
  esac
done

# ---------------------------------------------------------------------------
# Logging helpers
# ---------------------------------------------------------------------------
START_TS=$(date +%s)
log()  { printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*"; }
die()  { printf '[%s] ERROR: %s\n' "$(date +%H:%M:%S)" "$*" >&2; exit 1; }
stage_start() { STAGE_NAME="$1"; STAGE_TS=$(date +%s); log "=== $STAGE_NAME ==="; }
stage_end()   { log "--- $STAGE_NAME done in $(( $(date +%s) - STAGE_TS ))s"; }

# ---------------------------------------------------------------------------
# Preflight
# ---------------------------------------------------------------------------
log "preflight"

command -v csv_maker      >/dev/null || die "csv_maker not on PATH (build the image, or add the build dir to PATH)"
command -v csv_maker_mc   >/dev/null || die "csv_maker_mc not on PATH"
command -v python3        >/dev/null || die "python3 not on PATH"

# Prefer whichever python has calc_xsections importable.
PY=""
for cand in "${PYTHON:-}" python3 /opt/venv/bin/python; do
  [[ -z "$cand" ]] && continue
  if command -v "$cand" >/dev/null 2>&1 && "$cand" -c "import calc_xsections" >/dev/null 2>&1; then
    PY="$cand"; break
  fi
done
[[ -n "$PY" ]] || die "no python with calc_xseconds importable; set PYTHON=/path/to/python"
log "python: $PY"

if ! "$PY" -c "import calc_xsections as c" >/dev/null 2>&1; then
  die "importing calc_xseconds raised; run it directly to see why"
fi
if [[ "$("$PY" -c 'import calc_xsections as c; print(1 if c.maid else 0)')" != "1" ]]; then
  die "maid_interface is not loaded, so the phi bin-centering correction would
     silently be set to 1.0 (a ~1.8x effect on some phi bins).
     Build it, or unset the requirement by accepting that approximation."
fi
log "maid_interface: loaded"

# The Fortran resolves its parameter tables through $CLAS_PARMS/spp_tbl.
if [[ -z "${CLAS_PARMS:-}" ]]; then
  for cand in /opt/e1d "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/exe/maid_xsection"; do
    if [[ -d "$cand/spp_tbl" ]]; then export CLAS_PARMS="$cand"; break; fi
  done
fi
[[ -n "${CLAS_PARMS:-}" && -d "${CLAS_PARMS}/spp_tbl" ]] \
  || die "CLAS_PARMS is not set to a directory containing spp_tbl/"
log "CLAS_PARMS: $CLAS_PARMS"

[[ -n "$MOMCORR" ]] || die "--momcorr is required (the analysis exits 99 without MOM_CORR)"
[[ -d "$MOMCORR" ]] || die "MOM_CORR directory not found: $MOMCORR"
export MOM_CORR="$MOMCORR"
export NUM_THREADS="$THREADS"
log "MOM_CORR: $MOMCORR"
log "NUM_THREADS: $THREADS"

check_dir() {  # label, dir
  [[ -n "$2" ]] || return 0
  [[ -d "$2" ]] || die "$1 directory not found: $2"
  local n
  n=$(find "$2" -maxdepth 1 -type f | wc -l | tr -d ' ')
  [[ "$n" -gt 0 ]] || die "$1 directory has no files: $2"
  log "$1: $2 ($n files)"
}
check_dir "exp" "$EXP"
[[ "$ONLY_EXPERIMENTAL" == "1" ]] || check_dir "sim" "$SIM"

mkdir -p "$OUT"
OUT="$(cd "$OUT" && pwd)"
DATA_CSV="$OUT/data_e1d.csv"
MC_CSV="$OUT/mc_e1d.csv"
COARSE_CSV="$OUT/coarse_xsec.csv"

# Warn about disk: the MC CSV is ~4.7 MB per simulation file (987 files -> 4.6 GB)
# and dominates the run.
avail_gb=$(df -Pk "$OUT" | awk 'NR==2{printf "%.1f", $4/1048576}')
if [[ -n "$SIM" && "$CUT" == "0" ]]; then
  nsim=$(find "$SIM" -maxdepth 1 -type f | wc -l | tr -d ' ')
  est=$(awk -v n="$nsim" 'BEGIN{printf "%.1f", n*4.7/1024}')
  log "disk: ${avail_gb} GiB free in $OUT; MC CSV from $nsim sim files is roughly ${est} GiB"
  awk -v a="$avail_gb" -v e="$est" 'BEGIN{exit !(a < e*1.3)}' \
    && die "not enough free space: need about ${est} GiB plus headroom, have ${avail_gb} GiB"
fi

# ---------------------------------------------------------------------------
# Stage 1: experimental data -> CSV
# ---------------------------------------------------------------------------
if [[ -n "$EXP" ]]; then
  stage_start "csv_maker (experimental)"
  if [[ -f "$DATA_CSV" && "$FORCE" == "0" ]]; then
    log "exists, skipping (--force to redo): $DATA_CSV"
  else
    cutargs=()
    [[ "$CUT" != "0" ]] && cutargs=(-cut "$CUT")
    csv_maker "${cutargs[@]}" -e1d "$EXP" -o "$OUT/data"
    [[ -f "$DATA_CSV" ]] || die "expected $DATA_CSV to exist after csv_maker"
    log "rows: $(wc -l < "$DATA_CSV" | tr -d ' ')"
  fi
  stage_end
fi

# ---------------------------------------------------------------------------
# Stage 2: simulation -> CSV (both mc_rec and thrown rows in one file)
# ---------------------------------------------------------------------------
if [[ -n "$SIM" && "$ONLY_EXPERIMENTAL" == "0" ]]; then
  stage_start "csv_maker_mc (simulation)"
  if [[ -f "$MC_CSV" && "$FORCE" == "0" ]]; then
    log "exists, skipping (--force to redo): $MC_CSV"
  else
    cutargs=()
    [[ "$CUT" != "0" ]] && cutargs=(-cut "$CUT")
    csv_maker_mc "${cutargs[@]}" -e1d "$SIM" -o "$OUT/mc"
    [[ -f "$MC_CSV" ]] || die "expected $MC_CSV to exist after csv_maker_mc"
    log "rows: $(wc -l < "$MC_CSV" | tr -d ' ')"
    log "size: $(du -h "$MC_CSV" | cut -f1)"
  fi
  stage_end
fi

# ---------------------------------------------------------------------------
# Stage 3: phi-resolved cross sections
# ---------------------------------------------------------------------------
XSEC_DIR="$OUT/crossSections"
if [[ "$SKIP_XSEC" == "1" ]]; then
  log "skipping crossSectionsCSV.py (--skip-xsec)"
elif [[ -z "$MC_CSV" || -z "$DATA_CSV" ]]; then
  log "skipping crossSectionsCSV.py: needs both $DATA_CSV and $MC_CSV"
elif [[ -f "$XSEC_DIR/full_results_${PHIBINS}.csv" && "$FORCE" == "0" ]]; then
  log "exists, skipping (--force to redo): $XSEC_DIR/full_results_${PHIBINS}.csv"
else
  stage_start "crossSectionsCSV.py (phi-resolved, $PHIBINS phi bins)"
  mkdir -p "$OUT/crossSections"
  set +e
  "$PY" "$(dirname "${BASH_SOURCE[0]}")/../python/crossSectionsCSV.py" \
      --mc "$MC_CSV" --data "$DATA_CSV" --folder "$OUT" --phibins "$PHIBINS" \
      2>&1 | tee "$OUT/crossSections.log" | tail -20
  rc=${PIPESTATUS[0]}
  set -e
  stage_end
  if [[ "$rc" -ne 0 ]]; then
    die "crossSectionsCSV.py failed (rc=$rc); see $OUT/crossSections.log"
  fi
  log "NOTE: at the current MC statistics this binning is not statistically"
  log "      determined (median error/value ~300, no bin with error < value)."
  log "      Use the coarse cross section from the next stage for results."
fi

# ---------------------------------------------------------------------------
# Stage 4: coarse cross section + per-bin MC statistics
# ---------------------------------------------------------------------------
if [[ -z "$MC_CSV" || -z "$DATA_CSV" ]]; then
  die "coarse binning needs both $DATA_CSV and $MC_CSV"
fi
stage_start "coarse_xsec.py (coarse binning + MC statistics)"
if [[ -f "$COARSE_CSV" && "$FORCE" == "0" ]]; then
  log "exists, skipping (--force to redo): $COARSE_CSV"
else
  "$PY" "$(dirname "${BASH_SOURCE[0]}")/../analysis/coarse_xsec.py" \
      "$MC_CSV" "$DATA_CSV" "$COARSE_CSV" 2>&1 | tee "$OUT/coarse_xsec.log" | tail -15
  [[ -f "$COARSE_CSV" ]] || die "expected $COARSE_CSV to exist after coarse_xsec.py"
fi
stage_end

# ---------------------------------------------------------------------------
# Stage 5: plots
# ---------------------------------------------------------------------------
stage_start "plots"
PLOT_DIR="$OUT/plots"
mkdir -p "$PLOT_DIR"
if "$PY" "$(dirname "${BASH_SOURCE[0]}")/../analysis/plot_coarse_xsec.py" \
     "$COARSE_CSV" "$PLOT_DIR" 2>&1 | tail -5; then
  log "wrote coarse cross-section plot"
else
  die "plot_coarse_xsec.py failed"
fi
stage_end

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
log "============================================================"
log "chain complete in $(( $(date +%s) - START_TS ))s"
log "outputs in $OUT:"
for f in "$DATA_CSV" "$MC_CSV" "$COARSE_CSV" \
         "$XSEC_DIR/full_results_${PHIBINS}.csv" \
         "$PLOT_DIR/coarse_xsec_full_mc.csv" "$PLOT_DIR/xsec_coarse_W_full_mc.png"; do
  [[ -e "$f" ]] && printf '  %10s  %s\n' "$(du -h "$f" | cut -f1)" "$f"
done
log "============================================================"
