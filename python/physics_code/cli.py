"""CLI tools for CLAS12 analysis."""

from __future__ import annotations

import typer
from pathlib import Path
from typing import Optional
from rich.console import Console
from rich.table import Table

app = typer.Typer(name="physics-code", help="CLAS12 electron-scattering analysis framework")
console = Console()

EXPERIMENT_BEAM_ENERGIES = {
    "e1d": 4.81726,
    "e1f": 5.479,
    "e16": 5.76959,
}


@app.command()
def analyze(
    input_path: str = typer.Argument(..., help="Input ROOT file or directory containing ROOT files"),
    experiment: str = typer.Option("e1d", help="Experiment (e1d, e1f, e16)"),
    output: str = typer.Option("output.parquet", help="Output file"),
    beam_energy: Optional[float] = typer.Option(None, help="Beam energy in GeV (overrides experiment default)"),
    mc: bool = typer.Option(False, "--mc", help="Enable MC (thrown kinematics)"),
    num_threads: int = typer.Option(0, help="Number of worker threads (0 = auto-detect)"),
    batch_size: int = typer.Option(16, help="Number of files to process concurrently"),
    csv: bool = typer.Option(False, "--csv", help="Output as CSV instead of Parquet"),
):
    """Analyze ROOT files and produce a Parquet file with event data."""
    from physics_code._lib import process_files_to_parquet

    if beam_energy is None:
        beam_energy = EXPERIMENT_BEAM_ENERGIES.get(experiment.lower())
        if beam_energy is None:
            console.print(f"[red]Unknown experiment: {experiment}[/red]")
            raise typer.Exit(1)

    # Resolve input files
    input_path = Path(input_path)
    if input_path.is_dir():
        input_files = sorted(str(p) for p in input_path.glob("*.root"))
        if not input_files:
            console.print(f"[red]No .root files found in {input_path}[/red]")
            raise typer.Exit(1)
    elif input_path.is_file():
        input_files = [str(input_path)]
    else:
        console.print(f"[red]Path not found: {input_path}[/red]")
        raise typer.Exit(1)

    threads_info = f"threads={num_threads}" if num_threads > 0 else "threads=auto"
    output_format = "csv" if csv else "parquet"
    console.print(f"[green]Processing {len(input_files)} files with {experiment} cuts (E_beam = {beam_energy} GeV){' [MC]' if mc else ''} ({threads_info}, batch_size={batch_size})[/green]")

    n_events = process_files_to_parquet(input_files, experiment.lower(), beam_energy, output, mc, num_threads, batch_size, output_format)

    console.print(f"[green]Wrote {n_events} events to {output}[/green]")


@app.command()
def info(
    input_file: str = typer.Argument(..., help="Input ROOT file"),
):
    """Show information about a ROOT file."""
    from physics_code._lib import read_root_file_py

    n_events = read_root_file_py(input_file)
    console.print(f"[green]File: {input_file}[/green]")
    console.print(f"  Events: {n_events}")


@app.command()
def plot(
    input_file: str = typer.Argument(..., help="Input Parquet file"),
    output_dir: str = typer.Option("plots", help="Output directory for plots"),
    hist: str = typer.Option("w,q2", help="Comma-separated list of histograms to produce"),
    bins: int = typer.Option(100, help="Number of bins"),
):
    """Plot histograms from a Parquet file."""
    import polars as pl
    import matplotlib.pyplot as plt
    from pathlib import Path

    df = pl.read_parquet(input_file)
    out_path = Path(output_dir)
    out_path.mkdir(exist_ok=True)

    hist_names = [h.strip() for h in hist.split(",")]

    for h in hist_names:
        if h not in df.columns:
            console.print(f"[yellow]Column '{h}' not found in data, skipping[/yellow]")
            continue

        data = df[h].to_numpy()
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.hist(data, bins=bins, alpha=0.7, edgecolor="black")
        ax.set_xlabel(h)
        ax.set_ylabel("Counts")
        ax.set_title(f"{h} Distribution")

        outfile = out_path / f"{h}.png"
        fig.savefig(outfile, dpi=150, bbox_inches="tight")
        plt.close(fig)
        console.print(f"[green]Saved {outfile}[/green]")


@app.command()
def plot_sector(
    input_file: str = typer.Argument(..., help="Input Parquet file"),
    output_dir: str = typer.Option("plots", help="Output directory for plots"),
    variable: str = typer.Option("w", help="Variable to plot by sector"),
    bins: int = typer.Option(100, help="Number of bins"),
):
    """Plot histograms split by sector."""
    import polars as pl
    import matplotlib.pyplot as plt
    from pathlib import Path

    df = pl.read_parquet(input_file)
    out_path = Path(output_dir)
    out_path.mkdir(exist_ok=True)

    if variable not in df.columns:
        console.print(f"[red]Column '{variable}' not found[/red]")
        raise typer.Exit(1)

    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()

    for sector in range(1, 7):
        ax = axes[sector - 1]
        sector_data = df.filter(pl.col("sector") == sector)[variable].to_numpy()
        if len(sector_data) > 0:
            ax.hist(sector_data, bins=bins, alpha=0.7, edgecolor="black")
        ax.set_title(f"Sector {sector}")
        ax.set_xlabel(variable)
        ax.set_ylabel("Counts")

    fig.suptitle(f"{variable} by Sector", fontsize=14)
    fig.tight_layout()

    outfile = out_path / f"{variable}_by_sector.png"
    fig.savefig(outfile, dpi=150, bbox_inches="tight")
    plt.close(fig)
    console.print(f"[green]Saved {outfile}[/green]")


if __name__ == "__main__":
    app()
