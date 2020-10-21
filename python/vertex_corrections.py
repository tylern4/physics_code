#!/usr/bin/env python
# coding: utf-8
import matplotlib  # noqa
matplotlib.use("agg")  # noqa

import argparse
import functools
import json
import operator
from typing import Dict

import boost_histogram as bh
import h10
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from tqdm import tqdm


def good_event(event) -> bool:
    """good_event

    ### Args:
        event (h10event): Event from h10 reader

    ### Returns:
        bool: Returns whether the event has a good electron
    """
    if event.dc_sect.size == 0 or event.dc_vz.size == 0:
        return False
    if event.p.size == 0:
        return False

    return True


def plot_corretions(vz: Dict[int, bh.Histogram], *, folder="plots", **kwargs) -> Dict[int, float]:
    """Plots vertex uncorrected and vertex corrected hists and retruns the correction

    Args:
        vz (Dict[int, bh.Histogram]): Dict of sector to vertex histograms

    Returns:
        Dict[int, float]: Dict of sector to correction factor for that sector
    """
    # Output dictionary
    corrections = dict()
    # Create fig for plot
    fig, ax = plt.subplots(figsize=[12, 9])

    # For each sector in 1 to 6
    for sec in range(1, 7):
        # Get results of np.histogram from boost histogram
        y, x = vz[sec].to_numpy()
        # Normalize the y values
        y = y / np.max(y)
        # Get bin centers for each point
        x = (x[1:] + x[:-1]) / 2.0

        # Find the point where the leading edge of the peak comes in
        # Look for a place where y is above 0.2 and then get the first ([0][0]) value of that
        x_val = np.where(y >= 0.2)[0][0]
        # Set results in dict
        corrections[sec] = np.round(x[x_val], decimals=3)

        # Plot the corrected histogram for each sector
        ax.errorbar(x-corrections[sec], y, yerr=stats.sem(y), marker=".",
                    linestyle="", label=f'Sec {sec}: {x[x_val]}')

    fig.legend()
    fig.savefig("plots/corrected.png")

    return corrections


def plot_2d(vz_vs_phi: Dict[int, bh.Histogram]) -> None:
    for sec in range(1, 7):
        # Compute the areas of each bin
        areas = functools.reduce(operator.mul, vz_vs_phi[sec].axes.widths)

        # Compute the density
        density = vz_vs_phi[sec].view() / vz_vs_phi[sec].sum() / areas

        # Make the plot
        fig, ax = plt.subplots()
        mesh = ax.pcolormesh(*vz_vs_phi[sec].axes.edges.T, density.T)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Vertex Corrections")
    parser.add_argument(
        "--data",
        dest="data_file_path",
        type=str,
        help="Data runs location",
        required=True,
    )
    parser.add_argument(
        "--folder",
        dest="out_folder",
        type=str,
        help="Folder for plots",
        required=False,
        default="plots",
    )

    args = parser.parse_args()

    if args.data_file_path[-1] == "/":
        args.data_file_path = args.data_file_path[:-1]

    root_reader = h10.h10_data()
    root_reader.add(f"{args.data_file_path}/*.root")

    # Create dict to hold histograms
    vz = dict()
    vz_vs_phi = dict()

    # Place histograms into dictionary for each sector
    for sec in range(1, 7):
        vz[sec] = bh.Histogram(bh.axis.Regular(100, -6, 10))
        vz_vs_phi[sec] = bh.Histogram(bh.axis.Regular(100, -6, 10),
                                      bh.axis.Regular(100, -np.pi, np.pi))

    # Open a progress bar
    with tqdm(total=root_reader.num_entries) as pbar:
        # Loop over events in the data
        for event in root_reader:
            pbar.update(1)
            # Check if there is a good event and if not continue
            if not good_event(event):
                continue

            # Get the sector
            sector = event.dc_sect[0]

            # Fill sector for vz and vz_vs_phi
            vz[sector].fill(event.dc_vz[0])
            phi = np.arctan2(event.cx[0], event.cy[0])
            vz_vs_phi[sector].fill(event.dc_vz[0], phi)

    corrections = plot_corretions(vz=vz)
    with open('corrections.json', 'w') as outFile:
        json.dump(corrections, outFile)
