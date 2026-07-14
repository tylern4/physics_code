"""Plotting utilities for CLAS12 analysis."""

from pathlib import Path
from typing import Optional

import numpy as np
import polars as pl


def load_data(filename: str) -> pl.DataFrame:
    """Load analysis results from a Parquet file."""
    return pl.read_parquet(filename)


def histogram(
    data: np.ndarray,
    bins: int = 100,
    range: Optional[tuple[float, float]] = None,
    xlabel: str = "",
    ylabel: str = "Counts",
    title: str = "",
    output: Optional[str] = None,
):
    """Create a histogram and optionally save to file."""
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.hist(data, bins=bins, range=range, alpha=0.7, edgecolor="black")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)

    if output:
        fig.savefig(output, dpi=150, bbox_inches="tight")
    plt.close(fig)

    return fig


def plot_w_q2(
    df: pl.DataFrame,
    w_bins: int = 100,
    q2_bins: int = 100,
    output: Optional[str] = None,
):
    """Plot W vs Q2 2D histogram."""
    import matplotlib.pyplot as plt

    w = df["w"].to_numpy()
    q2 = df["q2"].to_numpy()

    fig, ax = plt.subplots(figsize=(8, 6))
    h = ax.hist2d(w, q2, bins=[w_bins, q2_bins], cmap="viridis")
    ax.set_xlabel("W (GeV)")
    ax.set_ylabel("Q² (GeV²)")
    ax.set_title("W vs Q²")
    plt.colorbar(h[3], ax=ax, label="Counts")

    if output:
        fig.savefig(output, dpi=150, bbox_inches="tight")
    plt.close(fig)

    return fig


def plot_missing_mass(
    df: pl.DataFrame,
    sector: Optional[int] = None,
    bins: int = 100,
    output: Optional[str] = None,
):
    """Plot missing mass distribution."""
    import matplotlib.pyplot as plt

    if sector is not None:
        data = df.filter(pl.col("sector") == sector)["mm2"].to_numpy()
        title = f"Missing Mass² (Sector {sector})"
    else:
        data = df["mm2"].to_numpy()
        title = "Missing Mass²"

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.hist(data, bins=bins, alpha=0.7, edgecolor="black")
    ax.set_xlabel("MM² (GeV²)")
    ax.set_ylabel("Counts")
    ax.set_title(title)

    if output:
        fig.savefig(output, dpi=150, bbox_inches="tight")
    plt.close(fig)

    return fig


def plot_theta_star(
    df: pl.DataFrame,
    bins: int = 100,
    output: Optional[str] = None,
):
    """Plot theta_star distribution."""
    import matplotlib.pyplot as plt

    data = df["theta_star"].to_numpy()

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.hist(data, bins=bins, alpha=0.7, edgecolor="black")
    ax.set_xlabel("θ* (rad)")
    ax.set_ylabel("Counts")
    ax.set_title("θ* Distribution")

    if output:
        fig.savefig(output, dpi=150, bbox_inches="tight")
    plt.close(fig)

    return fig
