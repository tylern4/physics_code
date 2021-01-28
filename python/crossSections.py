#!/usr/bin/env python

import matplotlib  # noqa

matplotlib.use("agg")  # noqa
import warnings  # noqa

warnings.filterwarnings("ignore")  # noqa

import argparse
import inspect
import os
import time
import matplotlib.pyplot as plt
from tqdm import tqdm
from calc_xsections import *


def main(rec, mc_rec, mc_thrown, binning, out_folder="plots", bins=24, overlap=None):
    if not os.path.exists(f'{out_folder}/crossSections'):
        os.makedirs(f'{out_folder}/crossSections')
    # Make a set of values from 0 to 2Pi for plotting
    xs = np.linspace(0, 2 * np.pi, 250)
    for w in tqdm(binning["wbins"]):
        for q2 in binning["q2bins"]:
            for theta in binning["thetabins"]:
                fig = plt.figure(
                    figsize=(12, 9), constrained_layout=True)
                gs = fig.add_gridspec(2, 1, height_ratios=[2, 1])
                ax1 = fig.add_subplot(gs[0])
                ax2 = fig.add_subplot(gs[1], sharex=ax1)
                if overlap is not None:
                    df = pd.read_csv(overlap, index_col=0)
                    old_data = df[(df.W_min == w.left) & (
                        df.Q2_min == q2.left) & (df.cos_t == theta.left)]
                    ebar = ax1.errorbar(old_data.phi, old_data.y, yerr=old_data.yerr,
                                        marker='*', linestyle="",
                                        zorder=1, label=f"E-99-107",
                                        markersize=10, alpha=0.4, c='r')

                plot_maid_model(ax1, w, q2, theta, xs)

                # Cut data/mc for the w/q2/theta bin we're in
                data = rec[make_cuts(rec, w, q2, theta)].copy()
                data_mc = mc_rec[make_cuts(mc_rec, w, q2, theta)].copy()
                thrown = mc_thrown[make_cuts(mc_thrown, w, q2, theta)].copy()
                num_good = np.sum(data.cut_fid)
                if num_good < 6:
                    continue
                elif num_good <= 24:
                    bins = 10
                else:
                    bins = 24

                cut_fids = {
                    "Fiducial Cuts": 0,
                    # "All Data": 2,
                    # "Fid cuts False": 1,
                }
                for name, cuts in cut_fids.items():
                    if cuts == 0:
                        # Histogram the data for plotting
                        marker = 'o'
                        _data_y, _x = hist_data(
                            data[data.cut_fid], density=False, bins=bins)
                        _mc_rec_y, _ = hist_data(
                            data_mc[data_mc.cut_fid], density=False, bins=bins)
                    elif cuts == 1:
                        # Histogram the data for plotting
                        marker = "^"
                        _data_y, _x = hist_data(
                            data[~data.cut_fid], density=False, bins=bins)
                        _mc_rec_y, _ = hist_data(
                            data_mc[~data_mc.cut_fid], density=False, bins=bins)
                    else:
                        marker = 'd'
                        _data_y, _x = hist_data(
                            data, density=False, bins=bins)
                        _mc_rec_y, _ = hist_data(
                            data_mc, density=False, bins=bins)

                    _thrown_y, _ = hist_data(thrown, density=False, bins=bins)

                    # Remove points with 0 data count
                    cut = ~(_data_y == 0)
                    x = _x[cut]
                    mc_rec_y = _mc_rec_y[cut]
                    thrown_y = _thrown_y[cut]
                    data_y = _data_y[cut]

                    # Remove points with 0 mc_rec count and put 1???
                    mc_rec_y = np.where(mc_rec_y == 0, 1, mc_rec_y)
                    thrown_y = np.where(thrown_y == 0, 1, thrown_y)
                    acceptance = mc_rec_y / thrown_y

                    # Plot intergrated yeils to compare with/without fid cuts
                    try:
                        # ax2.errorbar(x, 1/acceptance, marker=marker,
                        #              linestyle="", zorder=1,
                        #              markersize=10, label=f"Counts: {name}", alpha=0.4)
                        ax2.errorbar(x, data_y, marker=marker,
                                     linestyle="", zorder=1,
                                     markersize=10, label=f"Counts: {name}", alpha=0.4)
                        ax2.legend(loc='upper right')
                    except ValueError:
                        pass

                    # Get bin widths
                    delta_W = (w.right-w.left)
                    delta_Q2 = (q2.right-q2.left)
                    delta_Theta = np.abs(theta.right-theta.left)
                    __phis = np.linspace(0, 2 * np.pi, bins)
                    delta_phi = __phis[1] - __phis[0]
                    kin_bin_width = delta_W * delta_Q2 * delta_Theta * delta_phi

                    # Normalize with bin widths
                    try:
                        data_y = data_y/kin_bin_width
                    except ValueError:
                        continue

                    # Calculate acceptance and correct data
                    flux = virtual_photon_flux(w.left, q2.left) * luminosity()

                    y = data_y / acceptance / flux

                    error_bar = get_error_bars(y, mc_rec_y, thrown_y)

                    ebar = ax1.errorbar(x, y, yerr=error_bar,
                                        marker=marker, linestyle="",
                                        zorder=1, label=f"{name}", markersize=10, alpha=0.4)
                    out = fit_model(ax1, model_new, x, y, xs,
                                    ebar[0].get_color(), name)
                    ax1.legend(loc='upper right')

                    if cuts == 0:
                        try:
                            top = np.max(y)*1.5
                        except ValueError:
                            top = np.nan
                        if np.isnan(top) or np.isinf(top):
                            ax1.set_ylim(bottom=0, top=1.0)
                        else:
                            ax1.set_ylim(bottom=0, top=top)

                ax1.set_title(f"$W$ : {w} , $Q^2$ : {q2} $\\theta$ : {theta}")
                fig.savefig(f"{out_folder}/crossSections/w_{w.left}_q2_{q2.left}_theta_{theta.left}.png",
                            bbox_inches='tight')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Make Cross Sections")
    parser.add_argument("--mc", dest="mc_data_file_path",
                        type=str, help="MC csv file", required=True)
    parser.add_argument("--data", dest="rec_data_file_path",
                        type=str, help="Data csv file", required=True)
    parser.add_argument("--folder", dest="out_folder", type=str,
                        help="Folder for plots", required=False, default="plots")
    parser.add_argument("--overlap", dest="overlap", type=str, help="Location of overlap data csv", required=False,
                        default=None)
    args = parser.parse_args()

    # Start to main

    # Specifically put in bin edges
    # TODO ##################### BINS ######################
    if args.overlap is None:
        w_bins = np.array([1.1, 1.125, 1.15, 1.175, 1.2, 1.225, 1.25, 1.275, 1.3,
                           1.325, 1.35, 1.375, 1.4, 1.425, 1.45, 1.475, 1.5, 1.525,
                           1.55, 1.575, 1.6, 1.625, 1.65, 1.675, 1.7, 1.725, 1.75,
                           1.775, 1.8])
        q2_bins = np.array([1.2, 1.6, 2.0, 2.4, 3.5])
        theta_bins = np.array([-1.0, -0.8, -0.6, -0.4, -0.2,
                               0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2])
    else:
        # TODO ##################### BINS that overlap with KPark ######################
        w_bins = np.array([1.1, 1.12, 1.14, 1.16, 1.18, 1.2, 1.22, 1.24, 1.26, 1.28, 1.3,
                           1.32, 1.34, 1.36, 1.38, 1.4, 1.42, 1.44, 1.46, 1.48, 1.5, 1.52,
                           1.54, 1.56, 1.58, 1.6, 1.62, 1.64, 1.66, 1.68, 1.7, 1.72, 1.74,
                           1.76, 1.78, 1.8])
        q2_bins = np.array([1.1, 1.30, 1.56, 1.87, 2.23, 2.66])
        theta_bins = np.array([-1.0, -0.8, -0.6, -0.4, -0.2,
                               0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2])
    # TODO ##################### BINS ######################

    print("Start setup")
    start = time.time_ns()
    # Load mc file
    mc_rec, mc_thrown = read_csv(args.mc_data_file_path)
    # Load reconstructed file
    rec = read_csv(args.rec_data_file_path, True)

    # Cut for missing mass
    rec, mc_rec = cut_for_MM(rec, mc_rec)

    # Make bins in the dataframes from the bins above
    rec = prep_for_ana(rec, w_bins, q2_bins, theta_bins)
    mc_rec = prep_for_ana(mc_rec, w_bins, q2_bins, theta_bins)
    mc_thrown = prep_for_ana(mc_thrown, w_bins, q2_bins, theta_bins)

    # Create dict of sorted bins
    binning = dict()
    binning["wbins"] = pd.Index.sort_values(pd.unique(rec.w_bin))
    binning["q2bins"] = pd.Index.sort_values(pd.unique(rec.q2_bin))
    binning["thetabins"] = pd.Index.sort_values(pd.unique(rec.theta_bin))
    end = time.time_ns()
    print(f"Done setup: {(end-start)/1E9:0.2f}Sec")

    main(rec, mc_rec, mc_thrown, binning,
         out_folder=args.out_folder, overlap=args.overlap)
