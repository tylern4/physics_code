"""Plotting utilities for CLAS12 analysis.

Recreates all histograms from the C++ physics_code using parquet data.
Each function corresponds to one or more C++ histograms.

Data is in long format: one row per particle per event.
Event-level columns (w, q2, etc.) are repeated for each particle in the event.
Particle-level columns (part_p, part_beta, etc.) vary per row.
"""

from pathlib import Path
from typing import Optional

import numpy as np
import polars as pl
import matplotlib.pyplot as plt


def load_data(filename: str) -> pl.DataFrame:
    """Load analysis results from a Parquet file."""
    return pl.read_parquet(filename)


def _save(fig, output: Optional[str]):
    if output:
        fig.savefig(output, dpi=150, bbox_inches="tight")
    plt.close(fig)


def _get_event_level(df: pl.DataFrame) -> pl.DataFrame:
    """Get one row per event (deduplicate on event-level columns)."""
    return df.unique(subset=["w", "q2", "e_sector", "event_type", "e_prime", "mm"])


# ============================================================================
# 1. W and Q² Histograms
# ============================================================================

def plot_w(df: pl.DataFrame, output: Optional[str] = None, bins: int = 500):
    """W distribution (all events). C++: W_hist"""
    ev = _get_event_level(df)
    fig, ax = plt.subplots()
    ax.hist(ev["w"].to_numpy(), bins=bins, range=(0, 3.25), alpha=0.7, edgecolor="black")
    ax.set_xlabel("W (GeV)")
    ax.set_ylabel("Counts")
    ax.set_title("W")
    _save(fig, output)


def plot_q2(df: pl.DataFrame, output: Optional[str] = None, bins: int = 500):
    """Q² distribution (all events). C++: Q2_hist"""
    ev = _get_event_level(df)
    fig, ax = plt.subplots()
    ax.hist(ev["q2"].to_numpy(), bins=bins, range=(0, 5), alpha=0.7, edgecolor="black")
    ax.set_xlabel("Q² (GeV²)")
    ax.set_ylabel("Counts")
    ax.set_title("Q²")
    _save(fig, output)


def plot_w_q2(df: pl.DataFrame, output: Optional[str] = None, wbins: int = 500, q2bins: int = 500):
    """W vs Q² 2D histogram (all events). C++: WvsQ2_hist"""
    ev = _get_event_level(df)
    fig, ax = plt.subplots()
    h = ax.hist2d(ev["w"].to_numpy(), ev["q2"].to_numpy(),
                  bins=[wbins, q2bins], range=[[0, 3.25], [0, 5]], cmap="viridis")
    ax.set_xlabel("W (GeV)")
    ax.set_ylabel("Q² (GeV²)")
    ax.set_title("W vs Q²")
    plt.colorbar(h[3], ax=ax, label="Counts")
    _save(fig, output)


def plot_w_by_sector(df: pl.DataFrame, output: Optional[str] = None, bins: int = 500):
    """W per sector. C++: W_sec_{1-6}"""
    ev = _get_event_level(df)
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    for sec in range(1, 7):
        ax = axes[(sec - 1) // 3, (sec - 1) % 3]
        data = ev.filter(pl.col("e_sector") == sec)["w"].to_numpy()
        ax.hist(data, bins=bins, range=(0, 3.25), alpha=0.7, edgecolor="black")
        ax.set_title(f"Sector {sec}")
        ax.set_xlabel("W (GeV)")
        ax.set_ylabel("Counts")
    fig.suptitle("W by Sector")
    fig.tight_layout()
    _save(fig, output)


def plot_w_q2_by_sector(df: pl.DataFrame, output: Optional[str] = None, wbins: int = 200, q2bins: int = 200):
    """W vs Q² per sector. C++: W_vs_Q2_sec_{1-6}"""
    ev = _get_event_level(df)
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    for sec in range(1, 7):
        ax = axes[(sec - 1) // 3, (sec - 1) % 3]
        data = ev.filter(pl.col("e_sector") == sec)
        h = ax.hist2d(data["w"].to_numpy(), data["q2"].to_numpy(),
                      bins=[wbins, q2bins], range=[[0, 3.25], [0, 5]], cmap="viridis")
        ax.set_title(f"Sector {sec}")
        ax.set_xlabel("W (GeV)")
        ax.set_ylabel("Q² (GeV²)")
        plt.colorbar(h[3], ax=ax, label="Counts")
    fig.suptitle("W vs Q² by Sector")
    fig.tight_layout()
    _save(fig, output)


def plot_w_q2_proton(df: pl.DataFrame, output: Optional[str] = None, wbins: int = 500, q2bins: int = 500):
    """W vs Q² for proton-tagged events. C++: WvsQ2_proton"""
    ev = _get_event_level(df).filter(pl.col("event_type") == 22)
    fig, ax = plt.subplots()
    h = ax.hist2d(ev["w"].to_numpy(), ev["q2"].to_numpy(),
                  bins=[wbins, q2bins], range=[[0, 3.25], [0, 5]], cmap="viridis")
    ax.set_xlabel("W (GeV)")
    ax.set_ylabel("Q² (GeV²)")
    ax.set_title("W vs Q² (proton)")
    plt.colorbar(h[3], ax=ax, label="Counts")
    _save(fig, output)


def plot_w_q2_pion(df: pl.DataFrame, output: Optional[str] = None, wbins: int = 500, q2bins: int = 500):
    """W vs Q² for π⁺ events. C++: WvsQ2_pion"""
    ev = _get_event_level(df).filter(pl.col("event_type") == 0)
    fig, ax = plt.subplots()
    h = ax.hist2d(ev["w"].to_numpy(), ev["q2"].to_numpy(),
                  bins=[wbins, q2bins], range=[[0, 3.25], [0, 5]], cmap="viridis")
    ax.set_xlabel("W (GeV)")
    ax.set_ylabel("Q² (GeV²)")
    ax.set_title("W vs Q² (π⁺)")
    plt.colorbar(h[3], ax=ax, label="Counts")
    _save(fig, output)


def plot_w_q2_neutron_pip(df: pl.DataFrame, output: Optional[str] = None, wbins: int = 500, q2bins: int = 500):
    """W vs Q² for nπ⁺ events. C++: WvsQ2_NeutronPip"""
    ev = _get_event_level(df).filter(pl.col("event_type") == 10)
    fig, ax = plt.subplots()
    h = ax.hist2d(ev["w"].to_numpy(), ev["q2"].to_numpy(),
                  bins=[wbins, q2bins], range=[[0, 3.25], [0, 5]], cmap="viridis")
    ax.set_xlabel("W (GeV)")
    ax.set_ylabel("Q² (GeV²)")
    ax.set_title("W vs Q² (nπ⁺)")
    plt.colorbar(h[3], ax=ax, label="Counts")
    _save(fig, output)


def plot_w_q2_channel(df: pl.DataFrame, output: Optional[str] = None, wbins: int = 500, q2bins: int = 500):
    """W vs Q² for π⁺N channel. C++: WvsQ2_channel"""
    ev = _get_event_level(df).filter(pl.col("event_type").is_in([0, 10]))
    fig, ax = plt.subplots()
    h = ax.hist2d(ev["w"].to_numpy(), ev["q2"].to_numpy(),
                  bins=[wbins, q2bins], range=[[0, 3.25], [0, 5]], cmap="viridis")
    ax.set_xlabel("W (GeV)")
    ax.set_ylabel("Q² (GeV²)")
    ax.set_title("W vs Q² (π⁺N channel)")
    plt.colorbar(h[3], ax=ax, label="Counts")
    _save(fig, output)


def plot_w_q2_p_pi0(df: pl.DataFrame, output: Optional[str] = None, wbins: int = 250, q2bins: int = 250):
    """W vs Q² for pπ⁰ events. C++: WvsQ2_Ppi0"""
    ev = _get_event_level(df).filter(pl.col("event_type") == 222)
    fig, ax = plt.subplots()
    h = ax.hist2d(ev["w"].to_numpy(), ev["q2"].to_numpy(),
                  bins=[wbins, q2bins], range=[[0, 3.25], [0, 5]], cmap="viridis")
    ax.set_xlabel("W (GeV)")
    ax.set_ylabel("Q² (GeV²)")
    ax.set_title("W vs Q² (pπ⁰)")
    plt.colorbar(h[3], ax=ax, label="Counts")
    _save(fig, output)


def plot_w_q2_elastic(df: pl.DataFrame, output: Optional[str] = None, wbins: int = 500, q2bins: int = 500):
    """W vs Q² for elastic events. C++: WvsQ2_elastic"""
    ev = _get_event_level(df).filter(
        (pl.col("event_type") == 22) & (pl.col("w") < 2.0)
    )
    fig, ax = plt.subplots()
    h = ax.hist2d(ev["w"].to_numpy(), ev["q2"].to_numpy(),
                  bins=[wbins, q2bins], range=[[0, 2.0], [0, 5]], cmap="viridis")
    ax.set_xlabel("W (GeV)")
    ax.set_ylabel("Q² (GeV²)")
    ax.set_title("W vs Q² (elastic)")
    plt.colorbar(h[3], ax=ax, label="Counts")
    _save(fig, output)


def plot_w_q2_binned(df: pl.DataFrame, output: Optional[str] = None):
    """W vs Q² binned (24×4). C++: WvsQ2_binned"""
    ev = _get_event_level(df).filter(pl.col("event_type").is_in([0, 10]))
    fig, ax = plt.subplots()
    h = ax.hist2d(ev["w"].to_numpy(), ev["q2"].to_numpy(),
                  bins=[24, 4], range=[[1.3, 1.8], [1.0, 4.5]], cmap="viridis")
    ax.set_xlabel("W (GeV)")
    ax.set_ylabel("Q² (GeV²)")
    ax.set_title("W vs Q² (binned)")
    plt.colorbar(h[3], ax=ax, label="Counts")
    _save(fig, output)


def plot_e_prime(df: pl.DataFrame, output: Optional[str] = None, bins: int = 500):
    """Scattered electron energy. C++: E_prime_hist"""
    ev = _get_event_level(df).filter(pl.col("event_type").is_in([0, 10]))
    fig, ax = plt.subplots()
    ax.hist(ev["e_prime"].to_numpy(), bins=bins, range=(0, 5), alpha=0.7, edgecolor="black")
    ax.set_xlabel("E' (GeV)")
    ax.set_ylabel("Counts")
    ax.set_title("E' (channel)")
    _save(fig, output)


def plot_q2_vs_xb(df: pl.DataFrame, output: Optional[str] = None, bins: int = 500):
    """Q² vs x_B. C++: Q2_vs_xb"""
    ev = _get_event_level(df).filter(pl.col("event_type").is_in([0, 10]))
    fig, ax = plt.subplots()
    h = ax.hist2d(ev["xb"].to_numpy(), ev["q2"].to_numpy(),
                  bins=[bins, bins], range=[[0.1, 0.6], [1.0, 3.5]], cmap="viridis")
    ax.set_xlabel("x_B")
    ax.set_ylabel("Q² (GeV²)")
    ax.set_title("Q² vs x_B")
    plt.colorbar(h[3], ax=ax, label="Counts")
    _save(fig, output)


# ============================================================================
# 2. Missing Mass Histograms
# ============================================================================

def plot_missing_mass(df: pl.DataFrame, output: Optional[str] = None, bins: int = 250):
    """Missing mass. C++: Missing_Mass"""
    pip = df.filter(pl.col("part_is_pip"))
    fig, ax = plt.subplots()
    ax.hist(pip["mm"].to_numpy(), bins=bins, range=(0, 3), alpha=0.7, edgecolor="black")
    ax.set_xlabel("MM (GeV)")
    ax.set_ylabel("Counts")
    ax.set_title("Missing Mass (π⁺)")
    _save(fig, output)


def plot_missing_mass_small(df: pl.DataFrame, output: Optional[str] = None, bins: int = 1000):
    """Missing mass (narrow). C++: Missing_Mass_small"""
    pip = df.filter(pl.col("part_is_pip"))
    fig, ax = plt.subplots()
    ax.hist(pip["mm"].to_numpy(), bins=bins, range=(0.8, 1.8), alpha=0.7, edgecolor="black")
    ax.set_xlabel("MM (GeV)")
    ax.set_ylabel("Counts")
    ax.set_title("Missing Mass (narrow)")
    _save(fig, output)


def plot_missing_mass_square(df: pl.DataFrame, output: Optional[str] = None, bins: int = 250):
    """Missing mass squared. C++: Missing_Mass_square"""
    pip = df.filter(pl.col("part_is_pip"))
    fig, ax = plt.subplots()
    ax.hist(pip["mm2"].to_numpy(), bins=bins, range=(0.7, 1.5), alpha=0.7, edgecolor="black")
    ax.set_xlabel("MM² (GeV²)")
    ax.set_ylabel("Counts")
    ax.set_title("Missing Mass² (π⁺)")
    _save(fig, output)


def plot_missing_mass_by_sector(df: pl.DataFrame, output: Optional[str] = None, bins: int = 250):
    """Missing mass per sector. C++: Missing_Mass_small_{0-5}"""
    pip = df.filter(pl.col("part_is_pip"))
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    for sec in range(1, 7):
        ax = axes[(sec - 1) // 3, (sec - 1) % 3]
        data = pip.filter(pl.col("part_sector") == sec)["mm"].to_numpy()
        ax.hist(data, bins=bins, range=(0.8, 1.3), alpha=0.7, edgecolor="black")
        ax.set_title(f"Sector {sec}")
        ax.set_xlabel("MM (GeV)")
        ax.set_ylabel("Counts")
    fig.suptitle("Missing Mass by Sector")
    fig.tight_layout()
    _save(fig, output)


def plot_mm2_by_sector(df: pl.DataFrame, output: Optional[str] = None, bins: int = 250):
    """Missing mass squared per sector. C++: Missing_Mass_Sq_small_{0-5}"""
    pip = df.filter(pl.col("part_is_pip"))
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    for sec in range(1, 7):
        ax = axes[(sec - 1) // 3, (sec - 1) % 3]
        data = pip.filter(pl.col("part_sector") == sec)["mm2"].to_numpy()
        ax.hist(data, bins=bins, range=(0.7, 1.5), alpha=0.7, edgecolor="black")
        ax.set_title(f"Sector {sec}")
        ax.set_xlabel("MM² (GeV²)")
        ax.set_ylabel("Counts")
    fig.suptitle("Missing Mass² by Sector")
    fig.tight_layout()
    _save(fig, output)


# ============================================================================
# 3. Momentum vs Beta (PID)
# ============================================================================

def plot_mom_vs_beta(df: pl.DataFrame, output: Optional[str] = None, bins: int = 500):
    """Momentum vs beta (all particles). C++: MomVsBeta"""
    fig, ax = plt.subplots()
    h = ax.hist2d(df["part_p"].to_numpy(), df["part_beta"].to_numpy(),
                  bins=[bins, bins], range=[[0, 5], [0.1, 1.2]], cmap="viridis")
    ax.set_xlabel("P (GeV)")
    ax.set_ylabel("β")
    ax.set_title("Momentum vs β (all)")
    plt.colorbar(h[3], ax=ax, label="Counts")
    _save(fig, output)


def plot_mom_vs_beta_pos(df: pl.DataFrame, output: Optional[str] = None, bins: int = 500):
    """Momentum vs beta (positive). C++: MomVsBeta_pos"""
    pos = df.filter(pl.col("part_q") == 1)
    fig, ax = plt.subplots()
    h = ax.hist2d(pos["part_p"].to_numpy(), pos["part_beta"].to_numpy(),
                  bins=[bins, bins], range=[[0, 5], [0.1, 1.2]], cmap="viridis")
    ax.set_xlabel("P (GeV)")
    ax.set_ylabel("β")
    ax.set_title("Momentum vs β (+)")
    plt.colorbar(h[3], ax=ax, label="Counts")
    _save(fig, output)


def plot_mom_vs_beta_neg(df: pl.DataFrame, output: Optional[str] = None, bins: int = 500):
    """Momentum vs beta (negative). C++: MomVsBeta_neg"""
    neg = df.filter(pl.col("part_q") == -1)
    fig, ax = plt.subplots()
    h = ax.hist2d(neg["part_p"].to_numpy(), neg["part_beta"].to_numpy(),
                  bins=[bins, bins], range=[[0, 5], [0.1, 1.2]], cmap="viridis")
    ax.set_xlabel("P (GeV)")
    ax.set_ylabel("β")
    ax.set_title("Momentum vs β (-)")
    plt.colorbar(h[3], ax=ax, label="Counts")
    _save(fig, output)


def plot_mom_vs_beta_neutral(df: pl.DataFrame, output: Optional[str] = None, bins: int = 500):
    """Momentum vs beta (neutral). C++: MomVsBeta_Fill_neutral"""
    neut = df.filter(pl.col("part_q") == 0)
    fig, ax = plt.subplots()
    h = ax.hist2d(neut["part_p"].to_numpy(), neut["part_beta"].to_numpy(),
                  bins=[bins, bins], range=[[0, 5], [0.1, 1.2]], cmap="viridis")
    ax.set_xlabel("P (GeV)")
    ax.set_ylabel("β")
    ax.set_title("Momentum vs β (neutral)")
    plt.colorbar(h[3], ax=ax, label="Counts")
    _save(fig, output)


def plot_mom_vs_beta_proton(df: pl.DataFrame, output: Optional[str] = None, bins: int = 500):
    """Momentum vs beta (protons). C++: MomVsBeta_proton_ID"""
    prot = df.filter(pl.col("part_is_prot"))
    fig, ax = plt.subplots()
    h = ax.hist2d(prot["part_p"].to_numpy(), prot["part_beta"].to_numpy(),
                  bins=[bins, bins], range=[[0, 5], [0.1, 1.2]], cmap="viridis")
    ax.set_xlabel("P (GeV)")
    ax.set_ylabel("β")
    ax.set_title("Momentum vs β (proton)")
    plt.colorbar(h[3], ax=ax, label="Counts")
    _save(fig, output)


def plot_mom_vs_beta_pip(df: pl.DataFrame, output: Optional[str] = None, bins: int = 500):
    """Momentum vs beta (π⁺). C++: MomVsBeta_Pi_ID"""
    pip = df.filter(pl.col("part_is_pip"))
    fig, ax = plt.subplots()
    h = ax.hist2d(pip["part_p"].to_numpy(), pip["part_beta"].to_numpy(),
                  bins=[bins, bins], range=[[0, 5], [0.1, 1.2]], cmap="viridis")
    ax.set_xlabel("P (GeV)")
    ax.set_ylabel("β")
    ax.set_title("Momentum vs β (π⁺)")
    plt.colorbar(h[3], ax=ax, label="Counts")
    _save(fig, output)


def plot_momentum(df: pl.DataFrame, output: Optional[str] = None, bins: int = 500):
    """Momentum distribution. C++: Momentum"""
    fig, ax = plt.subplots()
    ax.hist(df["part_p"].to_numpy(), bins=bins, range=(0, 2.5), alpha=0.7, edgecolor="black")
    ax.set_xlabel("P (GeV)")
    ax.set_ylabel("Counts")
    ax.set_title("Momentum")
    _save(fig, output)


# ============================================================================
# 4. Delta-t (Time-of-Flight) Histograms
# ============================================================================

def plot_dt_proton(df: pl.DataFrame, output: Optional[str] = None, bins: int = 500):
    """P vs Δt (proton mass). C++: delta_t_mass_P"""
    pos = df.filter(pl.col("part_q") == 1)
    fig, ax = plt.subplots()
    h = ax.hist2d(pos["part_p"].to_numpy(), pos["part_delta_t_p"].to_numpy(),
                  bins=[bins, bins], range=[[0, 5], [-10, 10]], cmap="viridis")
    ax.set_xlabel("P (GeV)")
    ax.set_ylabel("Δt (ns)")
    ax.set_title("Δt (proton mass, q=+1)")
    plt.colorbar(h[3], ax=ax, label="Counts")
    _save(fig, output)


def plot_dt_proton_pid(df: pl.DataFrame, output: Optional[str] = None, bins: int = 500):
    """P vs Δt (proton, PID proton). C++: delta_t_mass_P_PID"""
    prot = df.filter(pl.col("part_is_prot"))
    fig, ax = plt.subplots()
    h = ax.hist2d(prot["part_p"].to_numpy(), prot["part_delta_t_p"].to_numpy(),
                  bins=[bins, bins], range=[[0, 5], [-10, 10]], cmap="viridis")
    ax.set_xlabel("P (GeV)")
    ax.set_ylabel("Δt (ns)")
    ax.set_title("Δt (proton, PID proton)")
    plt.colorbar(h[3], ax=ax, label="Counts")
    _save(fig, output)


def plot_dt_pip(df: pl.DataFrame, output: Optional[str] = None, bins: int = 500):
    """P vs Δt (π⁺ mass). C++: delta_t_mass_PIP"""
    pos = df.filter(pl.col("part_q") == 1)
    fig, ax = plt.subplots()
    h = ax.hist2d(pos["part_p"].to_numpy(), pos["part_delta_t_pip"].to_numpy(),
                  bins=[bins, bins], range=[[0, 5], [-10, 10]], cmap="viridis")
    ax.set_xlabel("P (GeV)")
    ax.set_ylabel("Δt (ns)")
    ax.set_title("Δt (π⁺ mass, q=+1)")
    plt.colorbar(h[3], ax=ax, label="Counts")
    _save(fig, output)


def plot_dt_pip_pid(df: pl.DataFrame, output: Optional[str] = None, bins: int = 500):
    """P vs Δt (π⁺, PID π⁺). C++: delta_t_mass_PIP_PID"""
    pip = df.filter(pl.col("part_is_pip"))
    fig, ax = plt.subplots()
    h = ax.hist2d(pip["part_p"].to_numpy(), pip["part_delta_t_pip"].to_numpy(),
                  bins=[bins, bins], range=[[0, 5], [-10, 10]], cmap="viridis")
    ax.set_xlabel("P (GeV)")
    ax.set_ylabel("Δt (ns)")
    ax.set_title("Δt (π⁺, PID π⁺)")
    plt.colorbar(h[3], ax=ax, label="Counts")
    _save(fig, output)


def plot_dt_pim(df: pl.DataFrame, output: Optional[str] = None, bins: int = 500):
    """P vs Δt (π⁻ mass). C++: delta_t_mass_PIM"""
    neg = df.filter(pl.col("part_q") == -1)
    fig, ax = plt.subplots()
    h = ax.hist2d(neg["part_p"].to_numpy(), neg["part_delta_t_pip"].to_numpy(),
                  bins=[bins, bins], range=[[0, 5], [-10, 10]], cmap="viridis")
    ax.set_xlabel("P (GeV)")
    ax.set_ylabel("Δt (ns)")
    ax.set_title("Δt (π⁻ mass, q=-1)")
    plt.colorbar(h[3], ax=ax, label="Counts")
    _save(fig, output)


def plot_dt_electron(df: pl.DataFrame, output: Optional[str] = None, bins: int = 500):
    """P vs Δt (electron mass). C++: delta_t_mass_electron"""
    neg = df.filter(pl.col("part_q") == -1)
    fig, ax = plt.subplots()
    h = ax.hist2d(neg["part_p"].to_numpy(), neg["part_delta_t_e"].to_numpy(),
                  bins=[bins, bins], range=[[0, 5], [-10, 10]], cmap="viridis")
    ax.set_xlabel("P (GeV)")
    ax.set_ylabel("Δt (ns)")
    ax.set_title("Δt (electron mass, q=-1)")
    plt.colorbar(h[3], ax=ax, label="Counts")
    _save(fig, output)


def plot_dt_kaon(df: pl.DataFrame, output: Optional[str] = None, bins: int = 500):
    """P vs Δt (K⁺ mass). C++: delta_t_mass_kp"""
    pos = df.filter(pl.col("part_q") == 1)
    fig, ax = plt.subplots()
    h = ax.hist2d(pos["part_p"].to_numpy(), pos["part_delta_t_k"].to_numpy(),
                  bins=[bins, bins], range=[[0, 5], [-10, 10]], cmap="viridis")
    ax.set_xlabel("P (GeV)")
    ax.set_ylabel("Δt (ns)")
    ax.set_title("Δt (K⁺ mass, q=+1)")
    plt.colorbar(h[3], ax=ax, label="Counts")
    _save(fig, output)


def plot_dt_proton_slices(df: pl.DataFrame, output: Optional[str] = None, bins: int = 500):
    """Δt (proton) in 20 momentum slices. C++: delta_t_p_{0-19}"""
    pos = df.filter(pl.col("part_q") == 1)
    fig, axes = plt.subplots(4, 5, figsize=(20, 16))
    for i in range(20):
        ax = axes[i // 5, i % 5]
        p_min = i * 0.25
        p_max = (i + 1) * 0.25
        data = pos.filter(
            (pl.col("part_p") >= p_min) & (pl.col("part_p") < p_max)
        )["part_delta_t_p"].to_numpy()
        ax.hist(data, bins=bins, range=(-10, 10), alpha=0.7, edgecolor="black")
        ax.set_title(f"{p_min:.1f}-{p_max:.1f} GeV")
        ax.set_xlabel("Δt (ns)")
    fig.suptitle("Δt (proton) by Momentum Slice")
    fig.tight_layout()
    _save(fig, output)


def plot_dt_pip_slices(df: pl.DataFrame, output: Optional[str] = None, bins: int = 500):
    """Δt (π⁺) in 20 momentum slices. C++: delta_t_pip_{0-19}"""
    pos = df.filter(pl.col("part_q") == 1)
    fig, axes = plt.subplots(4, 5, figsize=(20, 16))
    for i in range(20):
        ax = axes[i // 5, i % 5]
        p_min = i * 0.25
        p_max = (i + 1) * 0.25
        data = pos.filter(
            (pl.col("part_p") >= p_min) & (pl.col("part_p") < p_max)
        )["part_delta_t_pip"].to_numpy()
        ax.hist(data, bins=bins, range=(-10, 10), alpha=0.7, edgecolor="black")
        ax.set_title(f"{p_min:.1f}-{p_max:.1f} GeV")
        ax.set_xlabel("Δt (ns)")
    fig.suptitle("Δt (π⁺) by Momentum Slice")
    fig.tight_layout()
    _save(fig, output)


def plot_dt_electron_slices(df: pl.DataFrame, output: Optional[str] = None, bins: int = 500):
    """Δt (electron) in 20 momentum slices. C++: delta_t_electron_{0-19}"""
    neg = df.filter(pl.col("part_q") == -1)
    fig, axes = plt.subplots(4, 5, figsize=(20, 16))
    for i in range(20):
        ax = axes[i // 5, i % 5]
        p_min = i * 0.25
        p_max = (i + 1) * 0.25
        data = neg.filter(
            (pl.col("part_p") >= p_min) & (pl.col("part_p") < p_max)
        )["part_delta_t_e"].to_numpy()
        ax.hist(data, bins=bins, range=(-10, 10), alpha=0.7, edgecolor="black")
        ax.set_title(f"{p_min:.1f}-{p_max:.1f} GeV")
        ax.set_xlabel("Δt (ns)")
    fig.suptitle("Δt (electron) by Momentum Slice")
    fig.tight_layout()
    _save(fig, output)


# ============================================================================
# 5. Cherenkov Counter (CC) Histograms
# ============================================================================

def plot_cc_nphe(df: pl.DataFrame, output: Optional[str] = None, bins: int = 50):
    """Photoelectrons distribution. C++: CC_sec{s}_both"""
    elec = df.filter(pl.col("part_is_electron"))
    fig, ax = plt.subplots()
    ax.hist(elec["part_nphe"].to_numpy(), bins=bins, range=(0, 250), alpha=0.7, edgecolor="black")
    ax.set_xlabel("N_{phe}")
    ax.set_ylabel("Counts")
    ax.set_title("CC Photoelectrons (electron)")
    _save(fig, output)


def plot_cc_nphe_by_sector(df: pl.DataFrame, output: Optional[str] = None, bins: int = 50):
    """Photoelectrons per sector. C++: CC_sec{s}_both"""
    elec = df.filter(pl.col("part_is_electron"))
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    for sec in range(1, 7):
        ax = axes[(sec - 1) // 3, (sec - 1) % 3]
        data = elec.filter(pl.col("part_sector") == sec)["part_nphe"].to_numpy()
        ax.hist(data, bins=bins, range=(0, 250), alpha=0.7, edgecolor="black")
        ax.set_title(f"Sector {sec}")
        ax.set_xlabel("N_{phe}")
        ax.set_ylabel("Counts")
    fig.suptitle("CC Photoelectrons by Sector")
    fig.tight_layout()
    _save(fig, output)


def plot_theta_vs_cc_segment(df: pl.DataFrame, output: Optional[str] = None):
    """CC segment vs θ_CC. C++: Theta_CC"""
    elec = df.filter(pl.col("part_is_electron"))
    fig, ax = plt.subplots()
    h = ax.hist2d(elec["part_cc_segm"].to_numpy(), elec["part_theta"].to_numpy(),
                  bins=[20, 60], range=[[0, 20], [0, 60]], cmap="viridis")
    ax.set_xlabel("CC Segment")
    ax.set_ylabel("θ_CC")
    ax.set_title("CC Segment vs θ")
    plt.colorbar(h[3], ax=ax, label="Counts")
    _save(fig, output)


def plot_cc_fiducial_xy(df: pl.DataFrame, output: Optional[str] = None, bins: int = 200):
    """CC y vs x (fiducial). C++: fid_cher_xy_{1-6}"""
    elec = df.filter(pl.col("part_is_electron"))
    fig, ax = plt.subplots()
    # CC x/y are approximated from cc_r and cc_theta/phi
    # Use DC SC position as proxy (already in parquet)
    h = ax.hist2d(elec["part_dc_ysc"].to_numpy(), elec["part_dc_xsc"].to_numpy(),
                  bins=[bins, bins], range=[[-150, 150], [0, 300]], cmap="viridis")
    ax.set_xlabel("CC y")
    ax.set_ylabel("CC x")
    ax.set_title("CC Fiducial xy")
    plt.colorbar(h[3], ax=ax, label="Counts")
    _save(fig, output)


def plot_cc_fiducial_xy_by_sector(df: pl.DataFrame, output: Optional[str] = None, bins: int = 200):
    """CC y vs x per sector. C++: fid_cher_xy_{1-6}"""
    elec = df.filter(pl.col("part_is_electron"))
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    for sec in range(1, 7):
        ax = axes[(sec - 1) // 3, (sec - 1) % 3]
        data = elec.filter(pl.col("part_sector") == sec)
        h = ax.hist2d(data["part_dc_ysc"].to_numpy(), data["part_dc_xsc"].to_numpy(),
                      bins=[bins, bins], range=[[-150, 150], [0, 300]], cmap="viridis")
        ax.set_title(f"Sector {sec}")
        ax.set_xlabel("CC y")
        ax.set_ylabel("CC x")
        plt.colorbar(h[3], ax=ax, label="Counts")
    fig.suptitle("CC Fiducial xy by Sector")
    fig.tight_layout()
    _save(fig, output)


def plot_photon_pair_mass(df: pl.DataFrame, output: Optional[str] = None, bins: int = 250):
    """Invariant mass of all photon pairs from individual photon rows.
    C++: Mass_pi0 (from Fill_Mass_photons)"""
    photons = df.filter(pl.col("part_id") == 22)
    if len(photons) == 0:
        return

    pair_masses = []
    # Group photons by event (using unique event identifiers)
    event_keys = ["w", "q2", "e_sector", "e_prime"]
    for _, event_group in photons.group_by(event_keys):
        if len(event_group) < 2:
            continue
        rows = event_group.to_dicts()
        for i in range(len(rows)):
            for j in range(i + 1, len(rows)):
                p1, p2 = rows[i], rows[j]
                theta1 = np.radians(p1["part_theta"])
                phi1 = np.radians(p1["part_phi"])
                px1 = p1["part_p"] * np.sin(theta1) * np.cos(phi1)
                py1 = p1["part_p"] * np.sin(theta1) * np.sin(phi1)
                pz1 = p1["part_p"] * np.cos(theta1)
                e1 = np.sqrt(px1**2 + py1**2 + pz1**2)

                theta2 = np.radians(p2["part_theta"])
                phi2 = np.radians(p2["part_phi"])
                px2 = p2["part_p"] * np.sin(theta2) * np.cos(phi2)
                py2 = p2["part_p"] * np.sin(theta2) * np.sin(phi2)
                pz2 = p2["part_p"] * np.cos(theta2)
                e2 = np.sqrt(px2**2 + py2**2 + pz2**2)

                m2 = (e1 + e2)**2 - (px1+px2)**2 - (py1+py2)**2 - (pz1+pz2)**2
                if m2 > 0:
                    pair_masses.append(np.sqrt(m2))

    if not pair_masses:
        return

    fig, ax = plt.subplots()
    ax.hist(pair_masses, bins=bins, range=(0, 0.5), alpha=0.7, edgecolor="black")
    ax.set_xlabel("M(γγ) (GeV)")
    ax.set_ylabel("Counts")
    ax.set_title("Photon Pair Invariant Mass")
    _save(fig, output)


# ============================================================================
# 6. Electron Fiducial Histograms
# ============================================================================

def plot_electron_fid(df: pl.DataFrame, output: Optional[str] = None, bins: int = 500):
    """φ vs θ (electron, all). C++: electron_fid"""
    elec = df.filter(pl.col("part_is_electron") & (pl.col("part_theta") > 5))
    fig, ax = plt.subplots()
    h = ax.hist2d(elec["part_phi"].to_numpy(), elec["part_theta"].to_numpy(),
                  bins=[bins, bins], range=[[-180, 180], [0, 80]], cmap="viridis")
    ax.set_xlabel("φ (deg)")
    ax.set_ylabel("θ (deg)")
    ax.set_title("Electron Fiducial (all)")
    plt.colorbar(h[3], ax=ax, label="Counts")
    _save(fig, output)


def plot_electron_fid_cut(df: pl.DataFrame, output: Optional[str] = None, bins: int = 500):
    """φ vs θ (electron, after fid cuts). C++: electron_fid_cut"""
    elec = df.filter(pl.col("part_elec_fid"))
    fig, ax = plt.subplots()
    h = ax.hist2d(elec["part_phi"].to_numpy(), elec["part_theta"].to_numpy(),
                  bins=[bins, bins], range=[[-180, 180], [0, 80]], cmap="viridis")
    ax.set_xlabel("φ (deg)")
    ax.set_ylabel("θ (deg)")
    ax.set_title("Electron Fiducial (cut)")
    plt.colorbar(h[3], ax=ax, label="Counts")
    _save(fig, output)


def plot_electron_fid_anti(df: pl.DataFrame, output: Optional[str] = None, bins: int = 500):
    """φ vs θ (electron, anti-cut). C++: electron_fid_anti"""
    elec = df.filter(pl.col("part_is_electron") & ~pl.col("part_elec_fid"))
    fig, ax = plt.subplots()
    h = ax.hist2d(elec["part_phi"].to_numpy(), elec["part_theta"].to_numpy(),
                  bins=[bins, bins], range=[[-180, 180], [0, 80]], cmap="viridis")
    ax.set_xlabel("φ (deg)")
    ax.set_ylabel("θ (deg)")
    ax.set_title("Electron Fiducial (anti-cut)")
    plt.colorbar(h[3], ax=ax, label="Counts")
    _save(fig, output)


def plot_electron_fid_by_sector(df: pl.DataFrame, output: Optional[str] = None, bins: int = 200):
    """φ vs θ per sector. C++: electron_fid_sec{1-6}"""
    elec = df.filter(pl.col("part_is_electron") & (pl.col("part_theta") > 5))
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    for sec in range(1, 7):
        ax = axes[(sec - 1) // 3, (sec - 1) % 3]
        data = elec.filter(pl.col("part_sector") == sec)
        h = ax.hist2d(data["part_phi"].to_numpy(), data["part_theta"].to_numpy(),
                      bins=[bins, bins], range=[[-180, 180], [0, 80]], cmap="viridis")
        ax.set_title(f"Sector {sec}")
        ax.set_xlabel("φ (deg)")
        ax.set_ylabel("θ (deg)")
        plt.colorbar(h[3], ax=ax, label="Counts")
    fig.suptitle("Electron Fiducial by Sector")
    fig.tight_layout()
    _save(fig, output)


# ============================================================================
# 7. DC XY Position (Fiducial)
# ============================================================================

def plot_dc_xy(df: pl.DataFrame, output: Optional[str] = None, bins: int = 500):
    """DC y_sc vs x_sc (all). C++: fid_dc_xy"""
    elec = df.filter(pl.col("part_is_electron"))
    fig, ax = plt.subplots()
    h = ax.hist2d(elec["part_dc_ysc"].to_numpy(), elec["part_dc_xsc"].to_numpy(),
                  bins=[bins, bins], range=[[-150, 150], [0, 300]], cmap="viridis")
    ax.set_xlabel("y_sc")
    ax.set_ylabel("x_sc")
    ax.set_title("DC xy (electron)")
    plt.colorbar(h[3], ax=ax, label="Counts")
    _save(fig, output)


def plot_dc_xy_cut(df: pl.DataFrame, output: Optional[str] = None, bins: int = 500):
    """DC y_sc vs x_sc (after fid cut). C++: fid_dc_xy_cut"""
    elec = df.filter(pl.col("part_elec_fid"))
    fig, ax = plt.subplots()
    h = ax.hist2d(elec["part_dc_ysc"].to_numpy(), elec["part_dc_xsc"].to_numpy(),
                  bins=[bins, bins], range=[[-150, 150], [0, 300]], cmap="viridis")
    ax.set_xlabel("y_sc")
    ax.set_ylabel("x_sc")
    ax.set_title("DC xy (electron, cut)")
    plt.colorbar(h[3], ax=ax, label="Counts")
    _save(fig, output)


def plot_dc_xy_by_sector(df: pl.DataFrame, output: Optional[str] = None, bins: int = 200):
    """DC xy per sector. C++: fid_dc_xy_{1-6}"""
    elec = df.filter(pl.col("part_is_electron"))
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    for sec in range(1, 7):
        ax = axes[(sec - 1) // 3, (sec - 1) % 3]
        data = elec.filter(pl.col("part_sector") == sec)
        h = ax.hist2d(data["part_dc_ysc"].to_numpy(), data["part_dc_xsc"].to_numpy(),
                      bins=[bins, bins], range=[[-150, 150], [0, 300]], cmap="viridis")
        ax.set_title(f"Sector {sec}")
        ax.set_xlabel("y_sc")
        ax.set_ylabel("x_sc")
        plt.colorbar(h[3], ax=ax, label="Counts")
    fig.suptitle("DC xy by Sector")
    fig.tight_layout()
    _save(fig, output)


# ============================================================================
# 8. Hadron Fiducial Histograms
# ============================================================================

def plot_hadron_fid(df: pl.DataFrame, output: Optional[str] = None, bins: int = 500):
    """φ vs θ (hadrons). C++: hadron_fid"""
    had = df.filter(pl.col("part_hadron_fid"))
    fig, ax = plt.subplots()
    h = ax.hist2d(had["part_phi"].to_numpy(), had["part_theta"].to_numpy(),
                  bins=[bins, bins], range=[[-180, 180], [0, 80]], cmap="viridis")
    ax.set_xlabel("φ (deg)")
    ax.set_ylabel("θ (deg)")
    ax.set_title("Hadron Fiducial")
    plt.colorbar(h[3], ax=ax, label="Counts")
    _save(fig, output)


def plot_proton_fid(df: pl.DataFrame, output: Optional[str] = None, bins: int = 500):
    """φ vs θ (protons). C++: proton_fid"""
    prot = df.filter(pl.col("part_is_prot"))
    fig, ax = plt.subplots()
    h = ax.hist2d(prot["part_phi"].to_numpy(), prot["part_theta"].to_numpy(),
                  bins=[bins, bins], range=[[-180, 180], [0, 80]], cmap="viridis")
    ax.set_xlabel("φ (deg)")
    ax.set_ylabel("θ (deg)")
    ax.set_title("Proton Fiducial")
    plt.colorbar(h[3], ax=ax, label="Counts")
    _save(fig, output)


def plot_pip_fid(df: pl.DataFrame, output: Optional[str] = None, bins: int = 500):
    """φ vs θ (π⁺). C++: pip_fid"""
    pip = df.filter(pl.col("part_is_pip"))
    fig, ax = plt.subplots()
    h = ax.hist2d(pip["part_phi"].to_numpy(), pip["part_theta"].to_numpy(),
                  bins=[bins, bins], range=[[-180, 180], [0, 80]], cmap="viridis")
    ax.set_xlabel("φ (deg)")
    ax.set_ylabel("θ (deg)")
    ax.set_title("π⁺ Fiducial")
    plt.colorbar(h[3], ax=ax, label="Counts")
    _save(fig, output)


def plot_hadron_fid_by_sector(df: pl.DataFrame, output: Optional[str] = None, bins: int = 200):
    """Hadron fiducial per sector (proton=0, π⁺=1). C++: hadron_fid_sec{1-6}_{0,1,2}"""
    for pid_name, pid_label in [("part_is_prot", "proton"), ("part_is_pip", "π⁺")]:
        had = df.filter(pl.col(pid_name))
        fig, axes = plt.subplots(2, 3, figsize=(15, 10))
        for sec in range(1, 7):
            ax = axes[(sec - 1) // 3, (sec - 1) % 3]
            data = had.filter(pl.col("part_sector") == sec)
            h = ax.hist2d(data["part_phi"].to_numpy(), data["part_theta"].to_numpy(),
                          bins=[bins, bins], range=[[-180, 180], [0, 80]], cmap="viridis")
            ax.set_title(f"Sector {sec}")
            ax.set_xlabel("φ (deg)")
            ax.set_ylabel("θ (deg)")
            plt.colorbar(h[3], ax=ax, label="Counts")
        fig.suptitle(f"{pid_label} Fiducial by Sector")
        fig.tight_layout()
        _save(fig, output.replace(".png", f"_{pid_label}.png") if output else None)


# ============================================================================
# 9. EC (Electromagnetic Calorimeter) Histograms
# ============================================================================

def plot_ec_sampling_fraction(df: pl.DataFrame, output: Optional[str] = None, bins: int = 500):
    """P vs etot/P. C++: EC_sampling_fraction"""
    elec = df.filter(pl.col("part_is_electron") & (pl.col("part_p") > 0))
    sf = (elec["part_etot"] / elec["part_p"]).to_numpy()
    fig, ax = plt.subplots()
    h = ax.hist2d(elec["part_p"].to_numpy(), sf,
                  bins=[bins, bins], range=[[0, 5], [0, 1]], cmap="viridis")
    ax.set_xlabel("P (GeV)")
    ax.set_ylabel("EC sampling fraction")
    ax.set_title("EC Sampling Fraction")
    plt.colorbar(h[3], ax=ax, label="Counts")
    _save(fig, output)


def plot_ec_inner_vs_outer(df: pl.DataFrame, output: Optional[str] = None, bins: int = 500):
    """EC_inner vs EC_outer. C++: ECin_ECout"""
    elec = df.filter(pl.col("part_is_electron"))
    fig, ax = plt.subplots()
    h = ax.hist2d(elec["part_ec_ei"].to_numpy(), elec["part_ec_eo"].to_numpy(),
                  bins=[bins, bins], range=[[0, 0.5], [0, 0.5]], cmap="viridis")
    ax.set_xlabel("EC_inner")
    ax.set_ylabel("EC_outer")
    ax.set_title("EC Inner vs Outer")
    plt.colorbar(h[3], ax=ax, label="Counts")
    _save(fig, output)


def plot_ec_tot_energy(df: pl.DataFrame, output: Optional[str] = None, bins: int = 500):
    """EC total energy. C++: EC_tot_energy"""
    elec = df.filter(pl.col("part_is_electron"))
    fig, ax = plt.subplots()
    ax.hist(elec["part_etot"].to_numpy(), bins=bins, range=(0, 1.5), alpha=0.7, edgecolor="black")
    ax.set_xlabel("E_{tot} (GeV)")
    ax.set_ylabel("Counts")
    ax.set_title("EC Total Energy")
    _save(fig, output)


def plot_ec_etot_vs_p(df: pl.DataFrame, output: Optional[str] = None, bins: int = 500):
    """P vs EC total energy. C++: EC_etot_vs_P"""
    elec = df.filter(pl.col("part_is_electron"))
    fig, ax = plt.subplots()
    h = ax.hist2d(elec["part_p"].to_numpy(), elec["part_etot"].to_numpy(),
                  bins=[bins, bins], range=[[0, 5], [0, 1.5]], cmap="viridis")
    ax.set_xlabel("P (GeV)")
    ax.set_ylabel("E_{tot} (GeV)")
    ax.set_title("EC E_{tot} vs P")
    plt.colorbar(h[3], ax=ax, label="Counts")
    _save(fig, output)


def plot_ec_sf_by_momentum(df: pl.DataFrame, output: Optional[str] = None, bins: int = 500):
    """SF per momentum slice. C++: ec_{0-19}, ec_cut_{0-19}"""
    elec = df.filter(pl.col("part_is_electron") & (pl.col("part_p") > 0))
    fig, axes = plt.subplots(4, 5, figsize=(20, 16))
    for i in range(20):
        ax = axes[i // 5, i % 5]
        p_min = i * 0.25
        p_max = (i + 1) * 0.25
        data = elec.filter(
            (pl.col("part_p") >= p_min) & (pl.col("part_p") < p_max)
        )
        if len(data) > 0:
            sf = (data["part_etot"] / data["part_p"]).to_numpy()
            ax.hist(sf, bins=bins, range=(0, 1), alpha=0.7, edgecolor="black")
        ax.set_title(f"{p_min:.1f}-{p_max:.1f} GeV")
        ax.set_xlabel("SF")
    fig.suptitle("EC Sampling Fraction by Momentum Slice")
    fig.tight_layout()
    _save(fig, output)


# ============================================================================
# 10. Beam Position & Target Vertex Histograms
# ============================================================================

def plot_beam_position(df: pl.DataFrame, output: Optional[str] = None, bins: int = 500):
    """DC vx vs vy (electron). C++: Beam_Position"""
    elec = df.filter(pl.col("part_is_electron"))
    fig, ax = plt.subplots()
    h = ax.hist2d(elec["e_dc_vx"].to_numpy(), elec["e_dc_vy"].to_numpy(),
                  bins=[bins, bins], range=[[-0.5, 0.5], [-0.5, 0.5]], cmap="viridis")
    ax.set_xlabel("v_x (cm)")
    ax.set_ylabel("v_y (cm)")
    ax.set_title("Beam Position")
    plt.colorbar(h[3], ax=ax, label="Counts")
    _save(fig, output)


def plot_beam_position_x(df: pl.DataFrame, output: Optional[str] = None, bins: int = 500):
    """DC vx. C++: Beam_Position_X"""
    elec = df.filter(pl.col("part_is_electron"))
    fig, ax = plt.subplots()
    ax.hist(elec["e_dc_vx"].to_numpy(), bins=bins, range=(-0.5, 0.5), alpha=0.7, edgecolor="black")
    ax.set_xlabel("v_x (cm)")
    ax.set_ylabel("Counts")
    ax.set_title("Beam Position X")
    _save(fig, output)


def plot_beam_position_y(df: pl.DataFrame, output: Optional[str] = None, bins: int = 500):
    """DC vy. C++: Beam_Position_Y"""
    elec = df.filter(pl.col("part_is_electron"))
    fig, ax = plt.subplots()
    ax.hist(elec["e_dc_vy"].to_numpy(), bins=bins, range=(-0.5, 0.5), alpha=0.7, edgecolor="black")
    ax.set_xlabel("v_y (cm)")
    ax.set_ylabel("Counts")
    ax.set_title("Beam Position Y")
    _save(fig, output)


def plot_beam_position_z(df: pl.DataFrame, output: Optional[str] = None, bins: int = 5000):
    """DC vz. C++: Beam_Position_Z"""
    elec = df.filter(pl.col("part_is_electron"))
    fig, ax = plt.subplots()
    ax.hist(elec["e_dc_vz"].to_numpy(), bins=bins, range=(-10, 15), alpha=0.7, edgecolor="black")
    ax.set_xlabel("v_z (cm)")
    ax.set_ylabel("Counts")
    ax.set_title("Beam Position Z")
    _save(fig, output)


def plot_target_vertex(df: pl.DataFrame, output: Optional[str] = None, bins: int = 500):
    """Target vertex x vs y. C++: Target_vertex_xy"""
    valid = df.filter(
        (pl.col("part_vx").abs() < 6) & (pl.col("part_vy").abs() < 6)
    )
    fig, ax = plt.subplots()
    h = ax.hist2d(valid["part_vx"].to_numpy(), valid["part_vy"].to_numpy(),
                  bins=[bins, bins], range=[[-6, 6], [-6, 6]], cmap="viridis")
    ax.set_xlabel("v_x (cm)")
    ax.set_ylabel("v_y (cm)")
    ax.set_title("Target Vertex xy")
    plt.colorbar(h[3], ax=ax, label="Counts")
    _save(fig, output)


def plot_target_vertex_z(df: pl.DataFrame, output: Optional[str] = None, bins: int = 100):
    """Target vertex z. C++: Target_vertex_Z"""
    valid = df.filter(pl.col("part_vz").abs() < 6)
    fig, ax = plt.subplots()
    ax.hist(valid["part_vz"].to_numpy(), bins=bins, range=(-6, 6), alpha=0.7, edgecolor="black")
    ax.set_xlabel("v_z (cm)")
    ax.set_ylabel("Counts")
    ax.set_title("Target Vertex Z")
    _save(fig, output)


# ============================================================================
# 11. Angular Histograms
# ============================================================================

def plot_theta_vs_phi(df: pl.DataFrame, output: Optional[str] = None, bins: int = 500):
    """θ vs φ (all particles). C++: ThetaVsPhi_hist"""
    fig, ax = plt.subplots()
    h = ax.hist2d(df["part_theta"].to_numpy(), df["part_phi"].to_numpy(),
                  bins=[bins, bins], range=[[0, 180], [0, 360]], cmap="viridis")
    ax.set_xlabel("θ (deg)")
    ax.set_ylabel("φ (deg)")
    ax.set_title("θ vs φ (all)")
    plt.colorbar(h[3], ax=ax, label="Counts")
    _save(fig, output)


def plot_theta_vs_phi_channel(df: pl.DataFrame, output: Optional[str] = None, bins: int = 100):
    """θ vs φ (π⁺N channel). C++: ThetaVsPhi_channel"""
    ev = df.filter(pl.col("event_type").is_in([0, 10]))
    fig, ax = plt.subplots()
    h = ax.hist2d(ev["part_theta"].to_numpy(), ev["part_phi"].to_numpy(),
                  bins=[bins, bins], range=[[0, 180], [0, 360]], cmap="viridis")
    ax.set_xlabel("θ (deg)")
    ax.set_ylabel("φ (deg)")
    ax.set_title("θ vs φ (π⁺N channel)")
    plt.colorbar(h[3], ax=ax, label="Counts")
    _save(fig, output)


def plot_cos_theta_star_vs_phi_star(df: pl.DataFrame, output: Optional[str] = None, bins: int = 500):
    """cos(θ*) vs φ* (all). C++: CosThetaVsPhi_hist"""
    fig, ax = plt.subplots()
    cos_theta = np.cos(df["theta_star"].to_numpy())
    h = ax.hist2d(cos_theta, df["phi_star"].to_numpy(),
                  bins=[bins, bins], range=[[-1, 1], [0, 6.28]], cmap="viridis")
    ax.set_xlabel("cos(θ*)")
    ax.set_ylabel("φ* (rad)")
    ax.set_title("cos(θ*) vs φ*")
    plt.colorbar(h[3], ax=ax, label="Counts")
    _save(fig, output)


def plot_cos_theta_star_vs_phi_star_channel(df: pl.DataFrame, output: Optional[str] = None, bins: int = 100):
    """cos(θ*) vs φ* (π⁺N channel). C++: CosThetaVsPhi_channel"""
    ev = _get_event_level(df).filter(pl.col("event_type").is_in([0, 10]))
    fig, ax = plt.subplots()
    cos_theta = np.cos(ev["theta_star"].to_numpy())
    h = ax.hist2d(cos_theta, ev["phi_star"].to_numpy(),
                  bins=[bins, bins], range=[[-1, 1], [0, 6.28]], cmap="viridis")
    ax.set_xlabel("cos(θ*)")
    ax.set_ylabel("φ* (rad)")
    ax.set_title("cos(θ*) vs φ* (π⁺N channel)")
    plt.colorbar(h[3], ax=ax, label="Counts")
    _save(fig, output)


# ============================================================================
# 12. Theta vs Momentum per Sector
# ============================================================================

def plot_theta_vs_p_by_sector(df: pl.DataFrame, output: Optional[str] = None, bins: int = 200):
    """P vs θ per sector (electron and π⁺). C++: elec_theta_p_{1-6}, pip_theta_p_{1-6}"""
    for pid_name, pid_label, col_prefix in [
        ("part_is_electron", "electron", "e"),
        ("part_is_pip", "π⁺", "pip"),
    ]:
        data = df.filter(pl.col(pid_name))
        fig, axes = plt.subplots(2, 3, figsize=(15, 10))
        for sec in range(1, 7):
            ax = axes[(sec - 1) // 3, (sec - 1) % 3]
            sec_data = data.filter(pl.col("part_sector") == sec)
            h = ax.hist2d(sec_data["part_p"].to_numpy(), sec_data["part_theta"].to_numpy(),
                          bins=[bins, bins], range=[[0, 5], [0, 70]], cmap="viridis")
            ax.set_title(f"Sector {sec}")
            ax.set_xlabel("P (GeV)")
            ax.set_ylabel("θ (deg)")
            plt.colorbar(h[3], ax=ax, label="Counts")
        fig.suptitle(f"{pid_label}: P vs θ by Sector")
        fig.tight_layout()
        _save(fig, output.replace(".png", f"_{col_prefix}.png") if output else None)


def plot_theta_star_vs_p_by_sector(df: pl.DataFrame, output: Optional[str] = None, bins: int = 200):
    """P vs Θ* per sector. C++: elec_theta_star_p_{1-6}, pip_theta_star_p_{1-6}"""
    ev = _get_event_level(df).filter(pl.col("event_type").is_in([0, 10]))
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    for sec in range(1, 7):
        ax = axes[(sec - 1) // 3, (sec - 1) % 3]
        sec_data = ev.filter(pl.col("e_sector") == sec)
        h = ax.hist2d(sec_data["e_prime"].to_numpy(), sec_data["theta_star"].to_numpy(),
                      bins=[bins, bins], range=[[0, 5], [0, 3.14]], cmap="viridis")
        ax.set_title(f"Sector {sec}")
        ax.set_xlabel("E' (GeV)")
        ax.set_ylabel("θ* (rad)")
        plt.colorbar(h[3], ax=ax, label="Counts")
    fig.suptitle("E' vs θ* by Sector")
    fig.tight_layout()
    _save(fig, output)


# ============================================================================
# 13. Energy Histograms
# ============================================================================

def plot_energy_no_cuts(df: pl.DataFrame, output: Optional[str] = None, bins: int = 500):
    """E' (no cuts). C++: Energy_no_cuts"""
    ev = _get_event_level(df)
    fig, ax = plt.subplots()
    data = ev.filter(pl.col("e_prime") > 0.1)["e_prime"].to_numpy()
    ax.hist(data, bins=bins, range=(0, 5), alpha=0.7, edgecolor="black")
    ax.set_xlabel("E' (GeV)")
    ax.set_ylabel("Counts")
    ax.set_title("E' (no cuts)")
    _save(fig, output)


def plot_energy_channel(df: pl.DataFrame, output: Optional[str] = None, bins: int = 500):
    """E' (channel events). C++: Energy_channel_cuts"""
    ev = _get_event_level(df).filter(pl.col("event_type").is_in([0, 10]))
    fig, ax = plt.subplots()
    data = ev.filter(pl.col("e_prime") > 0.1)["e_prime"].to_numpy()
    ax.hist(data, bins=bins, range=(0, 5), alpha=0.7, edgecolor="black")
    ax.set_xlabel("E' (GeV)")
    ax.set_ylabel("Counts")
    ax.set_title("E' (channel)")
    _save(fig, output)


# ============================================================================
# 14. MC Histograms
# ============================================================================

def plot_w_q2_mc(df: pl.DataFrame, output: Optional[str] = None, wbins: int = 500, q2bins: int = 500):
    """W vs Q² (MC thrown). C++: WvsQ2_MC"""
    ev = _get_event_level(df).filter(pl.col("w_thrown").is_not_nan())
    fig, ax = plt.subplots()
    h = ax.hist2d(ev["w_thrown"].to_numpy(), ev["q2_thrown"].to_numpy(),
                  bins=[wbins, q2bins], range=[[0, 3.25], [0, 5]], cmap="viridis")
    ax.set_xlabel("W_thrown (GeV)")
    ax.set_ylabel("Q²_thrown (GeV²)")
    ax.set_title("W vs Q² (MC thrown)")
    plt.colorbar(h[3], ax=ax, label="Counts")
    _save(fig, output)


def plot_w_mc(df: pl.DataFrame, output: Optional[str] = None, bins: int = 500):
    """W_thrown. C++: W_MC"""
    ev = _get_event_level(df).filter(pl.col("w_thrown").is_not_nan())
    fig, ax = plt.subplots()
    ax.hist(ev["w_thrown"].to_numpy(), bins=bins, range=(0, 3.25), alpha=0.7, edgecolor="black")
    ax.set_xlabel("W_thrown (GeV)")
    ax.set_ylabel("Counts")
    ax.set_title("W (MC thrown)")
    _save(fig, output)


def plot_w_q2_mc_binned(df: pl.DataFrame, output: Optional[str] = None):
    """W vs Q² binned (MC). C++: WvsQ2_hist_binned_MC"""
    ev = _get_event_level(df).filter(pl.col("w_thrown").is_not_nan())
    fig, ax = plt.subplots()
    h = ax.hist2d(ev["w_thrown"].to_numpy(), ev["q2_thrown"].to_numpy(),
                  bins=[24, 4], range=[[1.3, 1.8], [1.0, 4.5]], cmap="viridis")
    ax.set_xlabel("W_thrown (GeV)")
    ax.set_ylabel("Q²_thrown (GeV²)")
    ax.set_title("W vs Q² binned (MC)")
    plt.colorbar(h[3], ax=ax, label="Counts")
    _save(fig, output)


# ============================================================================
# 15. Elastic / P π⁰ Specific
# ============================================================================

def plot_elastic_mm2(df: pl.DataFrame, output: Optional[str] = None, bins: int = 500):
    """MM² (elastic). C++: Elastic_MM"""
    prot = df.filter(pl.col("part_is_prot"))
    fig, ax = plt.subplots()
    ax.hist(prot["mm2"].to_numpy(), bins=bins, range=(-0.2, 0.2), alpha=0.7, edgecolor="black")
    ax.set_xlabel("MM² (GeV²)")
    ax.set_ylabel("Counts")
    ax.set_title("MM² (proton, elastic)")
    _save(fig, output)


def plot_pi0_mass(df: pl.DataFrame, output: Optional[str] = None, bins: int = 250):
    """π⁰ pair mass from 2 photons. C++: Mass_pi0"""
    ev = _get_event_level(df).filter(pl.col("num_photons") == 2)
    data = ev.filter(pl.col("pi0_mass").is_not_nan())["pi0_mass"].to_numpy()
    fig, ax = plt.subplots()
    ax.hist(data, bins=bins, range=(0, 0.5), alpha=0.7, edgecolor="black")
    ax.set_xlabel("M(γγ) (GeV)")
    ax.set_ylabel("Counts")
    ax.set_title("π⁰ Mass (2γ)")
    _save(fig, output)


def plot_pi0_mass2(df: pl.DataFrame, output: Optional[str] = None, bins: int = 250):
    """π⁰ pair mass² from 2 photons. C++: Mass_pi0_2"""
    ev = _get_event_level(df).filter(pl.col("num_photons") == 2)
    data = ev.filter(pl.col("pi0_mass2").is_not_nan())["pi0_mass2"].to_numpy()
    fig, ax = plt.subplots()
    ax.hist(data, bins=bins, range=(0, 0.05), alpha=0.7, edgecolor="black")
    ax.set_xlabel("M²(γγ) (GeV²)")
    ax.set_ylabel("Counts")
    ax.set_title("π⁰ Mass² (2γ)")
    _save(fig, output)


def plot_missing_mass_pi0(df: pl.DataFrame, output: Optional[str] = None, bins: int = 250):
    """MM for p events (π⁰ candidate). C++: Missing_Mass_pi0"""
    prot = df.filter(pl.col("part_is_prot") & (df["w"] < 2.0))
    fig, ax = plt.subplots()
    ax.hist(prot["mm"].to_numpy(), bins=bins, range=(-1, 1), alpha=0.7, edgecolor="black")
    ax.set_xlabel("MM (GeV)")
    ax.set_ylabel("Counts")
    ax.set_title("Missing Mass (p, W<2)")
    _save(fig, output)


def plot_missing_mass_pi0_2(df: pl.DataFrame, output: Optional[str] = None, bins: int = 250):
    """MM² for p events (π⁰ candidate). C++: Missing_Mass_pi0_2"""
    prot = df.filter(pl.col("part_is_prot") & (df["w"] < 2.0))
    fig, ax = plt.subplots()
    ax.hist(prot["mm2"].to_numpy(), bins=bins, range=(-1, 1), alpha=0.7, edgecolor="black")
    ax.set_xlabel("MM² (GeV²)")
    ax.set_ylabel("Counts")
    ax.set_title("Missing Mass² (p, W<2)")
    _save(fig, output)


# ============================================================================
# 16. Missing Mass (Other Channels)
# ============================================================================

def plot_missing_mass_two_pion(df: pl.DataFrame, output: Optional[str] = None, bins: int = 250):
    """Missing mass (2π). C++: Missing_Mass_2pi"""
    ev = _get_event_level(df).filter(pl.col("event_type") == 3333)
    fig, ax = plt.subplots()
    ax.hist(ev["mm"].to_numpy(), bins=bins, range=(0, 3), alpha=0.7, edgecolor="black")
    ax.set_xlabel("MM (GeV)")
    ax.set_ylabel("Counts")
    ax.set_title("Missing Mass (2π)")
    _save(fig, output)


def plot_mm2_two_pion(df: pl.DataFrame, output: Optional[str] = None, bins: int = 250):
    """Missing mass² (2π). C++: Missing_Mass_square_2pi"""
    ev = _get_event_level(df).filter(pl.col("event_type") == 3333)
    fig, ax = plt.subplots()
    ax.hist(ev["mm2"].to_numpy(), bins=bins, range=(0, 9), alpha=0.7, edgecolor="black")
    ax.set_xlabel("MM² (GeV²)")
    ax.set_ylabel("Counts")
    ax.set_title("Missing Mass² (2π)")
    _save(fig, output)


# ============================================================================
# 17. Phi Difference (Elastic)
# ============================================================================

def plot_phi_diff(df: pl.DataFrame, output: Optional[str] = None, bins: int = 500):
    """φ difference (e-p). C++: Elastic_phi"""
    prot = df.filter(pl.col("part_is_prot"))
    fig, ax = plt.subplots()
    data = prot.filter(
        (prot["mm2"].abs() < 0.005)
    )["part_phi"].to_numpy()
    ax.hist(data, bins=bins, range=(2, 4), alpha=0.7, edgecolor="black")
    ax.set_xlabel("Δφ (deg)")
    ax.set_ylabel("Counts")
    ax.set_title("φ difference (elastic)")
    _save(fig, output)


# ============================================================================
# 18. W and Q² Channel Histograms
# ============================================================================

def plot_w_channel(df: pl.DataFrame, output: Optional[str] = None, bins: int = 500):
    """W (π⁺N channel). C++: W_channel"""
    ev = _get_event_level(df).filter(pl.col("event_type").is_in([0, 10]))
    fig, ax = plt.subplots()
    ax.hist(ev["w"].to_numpy(), bins=bins, range=(0, 3.25), alpha=0.7, edgecolor="black")
    ax.set_xlabel("W (GeV)")
    ax.set_ylabel("Counts")
    ax.set_title("W (π⁺N channel)")
    _save(fig, output)


def plot_q2_channel(df: pl.DataFrame, output: Optional[str] = None, bins: int = 500):
    """Q² (π⁺N channel). C++: Q2_channel"""
    ev = _get_event_level(df).filter(pl.col("event_type").is_in([0, 10]))
    fig, ax = plt.subplots()
    ax.hist(ev["q2"].to_numpy(), bins=bins, range=(0, 5), alpha=0.7, edgecolor="black")
    ax.set_xlabel("Q² (GeV²)")
    ax.set_ylabel("Counts")
    ax.set_title("Q² (π⁺N channel)")
    _save(fig, output)


# ============================================================================
# Comprehensive Plot-All
# ============================================================================

ALL_PLOTS = {
    # W and Q²
    "w": plot_w,
    "q2": plot_q2,
    "w_q2": plot_w_q2,
    "w_by_sector": plot_w_by_sector,
    "w_q2_by_sector": plot_w_q2_by_sector,
    "w_q2_proton": plot_w_q2_proton,
    "w_q2_pion": plot_w_q2_pion,
    "w_q2_neutron_pip": plot_w_q2_neutron_pip,
    "w_q2_channel": plot_w_q2_channel,
    "w_q2_p_pi0": plot_w_q2_p_pi0,
    "w_q2_elastic": plot_w_q2_elastic,
    "w_q2_binned": plot_w_q2_binned,
    "e_prime": plot_e_prime,
    "q2_vs_xb": plot_q2_vs_xb,
    # Missing Mass
    "missing_mass": plot_missing_mass,
    "missing_mass_small": plot_missing_mass_small,
    "missing_mass_sq": plot_missing_mass_square,
    "missing_mass_by_sector": plot_missing_mass_by_sector,
    "mm2_by_sector": plot_mm2_by_sector,
    # PID
    "mom_vs_beta": plot_mom_vs_beta,
    "mom_vs_beta_pos": plot_mom_vs_beta_pos,
    "mom_vs_beta_neg": plot_mom_vs_beta_neg,
    "mom_vs_beta_neutral": plot_mom_vs_beta_neutral,
    "mom_vs_beta_proton": plot_mom_vs_beta_proton,
    "mom_vs_beta_pip": plot_mom_vs_beta_pip,
    "momentum": plot_momentum,
    # Delta-t
    "dt_proton": plot_dt_proton,
    "dt_proton_pid": plot_dt_proton_pid,
    "dt_pip": plot_dt_pip,
    "dt_pip_pid": plot_dt_pip_pid,
    "dt_pim": plot_dt_pim,
    "dt_electron": plot_dt_electron,
    "dt_kaon": plot_dt_kaon,
    "dt_proton_slices": plot_dt_proton_slices,
    "dt_pip_slices": plot_dt_pip_slices,
    "dt_electron_slices": plot_dt_electron_slices,
    # CC
    "cc_nphe": plot_cc_nphe,
    "cc_nphe_by_sector": plot_cc_nphe_by_sector,
    "theta_vs_cc_segment": plot_theta_vs_cc_segment,
    "cc_fiducial_xy": plot_cc_fiducial_xy,
    "cc_fiducial_xy_by_sector": plot_cc_fiducial_xy_by_sector,
    # Fiducial
    "electron_fid": plot_electron_fid,
    "electron_fid_cut": plot_electron_fid_cut,
    "electron_fid_anti": plot_electron_fid_anti,
    "electron_fid_by_sector": plot_electron_fid_by_sector,
    "dc_xy": plot_dc_xy,
    "dc_xy_cut": plot_dc_xy_cut,
    "dc_xy_by_sector": plot_dc_xy_by_sector,
    "hadron_fid": plot_hadron_fid,
    "proton_fid": plot_proton_fid,
    "pip_fid": plot_pip_fid,
    "hadron_fid_by_sector": plot_hadron_fid_by_sector,
    # EC
    "ec_sampling_fraction": plot_ec_sampling_fraction,
    "ec_inner_vs_outer": plot_ec_inner_vs_outer,
    "ec_tot_energy": plot_ec_tot_energy,
    "ec_etot_vs_p": plot_ec_etot_vs_p,
    "ec_sf_by_momentum": plot_ec_sf_by_momentum,
    # Beam position
    "beam_position": plot_beam_position,
    "beam_position_x": plot_beam_position_x,
    "beam_position_y": plot_beam_position_y,
    "beam_position_z": plot_beam_position_z,
    "target_vertex": plot_target_vertex,
    "target_vertex_z": plot_target_vertex_z,
    # Angular
    "theta_vs_phi": plot_theta_vs_phi,
    "theta_vs_phi_channel": plot_theta_vs_phi_channel,
    "cos_theta_star_vs_phi_star": plot_cos_theta_star_vs_phi_star,
    "cos_theta_star_vs_phi_star_channel": plot_cos_theta_star_vs_phi_star_channel,
    # Theta vs P
    "theta_vs_p_by_sector": plot_theta_vs_p_by_sector,
    "theta_star_vs_p_by_sector": plot_theta_star_vs_p_by_sector,
    # Energy
    "energy_no_cuts": plot_energy_no_cuts,
    "energy_channel": plot_energy_channel,
    # MC
    "w_q2_mc": plot_w_q2_mc,
    "w_mc": plot_w_mc,
    "w_q2_mc_binned": plot_w_q2_mc_binned,
    # Elastic
    "elastic_mm2": plot_elastic_mm2,
    "phi_diff": plot_phi_diff,
    "pi0_mass": plot_pi0_mass,
    "pi0_mass2": plot_pi0_mass2,
    "missing_mass_pi0": plot_missing_mass_pi0,
    "missing_mass_pi0_2": plot_missing_mass_pi0_2,
    "photon_pair_mass": plot_photon_pair_mass,
    # Two pion
    "missing_mass_2pi": plot_missing_mass_two_pion,
    "mm2_2pi": plot_mm2_two_pion,
    # Channel
    "w_channel": plot_w_channel,
    "q2_channel": plot_q2_channel,
}


def plot_all(
    df: pl.DataFrame,
    output_dir: str = "plots",
    plot_names: Optional[list[str]] = None,
):
    """Generate all histograms from the parquet data.

    Args:
        df: Loaded parquet DataFrame
        output_dir: Directory to save plots
        plot_names: Optional list of specific plot names to generate.
                    If None, generates all plots.
    """
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    plots_to_run = plot_names if plot_names else list(ALL_PLOTS.keys())

    for name in plots_to_run:
        if name not in ALL_PLOTS:
            print(f"Warning: unknown plot '{name}', skipping")
            continue
        try:
            fn = ALL_PLOTS[name]
            fn(df, output=str(out_path / f"{name}.png"))
            print(f"  Saved {name}.png")
        except Exception as e:
            print(f"  Error generating {name}: {e}")
