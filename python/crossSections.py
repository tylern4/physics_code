#!/usr/bin/env python
import matplotlib  # noqa
matplotlib.use('agg')  # noqa
import warnings  # noqa
warnings.filterwarnings("ignore")  # noqa

from calc_xsections import *
from tqdm import tqdm
import argparse
import inspect
import os
import time
import matplotlib.pyplot as plt

plt.rcParams.update({'mathtext.fontset': 'stix'})


@np.vectorize
def isclose(a, b, rel_tol=1e-4, abs_tol=0.0):
    return np.abs(a-b) <= np.maximum(rel_tol * np.maximum(np.abs(a), np.abs(b)), abs_tol)


def main(rec, mc_rec, mc_thrown, empty, binning, out_folder="plots", bins=12, overlap=None, radcorr=None):
    if not os.path.exists(f'{out_folder}/crossSections'):
        os.makedirs(f'{out_folder}/crossSections')
    # Make a set of values from 0 to 2Pi for plotting
    xs = np.linspace(0, 2 * np.pi, 250)
    if overlap is not None:
        overlap_df = pd.read_csv(overlap)
    if radcorr is not None:
        radcorr_df = pd.read_csv(radcorr)

    for w in tqdm(binning["wbins"]):
        for q2 in binning["q2bins"]:

            CosTfig = plt.figure(
                figsize=(12, 9), constrained_layout=True)
            ct_gs = CosTfig.add_gridspec(5, 2, hspace=0.1)
            _left_ax = CosTfig.add_subplot(ct_gs[0, 0])
            _right_ax = CosTfig.add_subplot(ct_gs[0, 1])
            ct_ax = {
                -1.0: _left_ax,
                -0.8: _right_ax,
                -0.6: CosTfig.add_subplot(ct_gs[1, 0], sharex=_left_ax),
                -0.4: CosTfig.add_subplot(ct_gs[1, 1], sharex=_right_ax),
                -0.2: CosTfig.add_subplot(ct_gs[2, 0], sharex=_left_ax),
                0.0: CosTfig.add_subplot(ct_gs[2, 1], sharex=_right_ax),
                0.2: CosTfig.add_subplot(ct_gs[3, 0], sharex=_left_ax),
                0.4: CosTfig.add_subplot(ct_gs[3, 1], sharex=_right_ax),
                0.6: CosTfig.add_subplot(ct_gs[4, 0], sharex=_left_ax),
                0.8: CosTfig.add_subplot(ct_gs[4, 1], sharex=_right_ax),
            }
            pass_plotting = False

            radcor_R = 1.0
            if radcorr is not None:
                cut = isclose(radcorr_df.w_left, w.left) & isclose(
                    radcorr_df.q2_left, q2.left)

                if(radcorr_df[cut].R.size == 0):
                    radcor_R = 1.0
                else:
                    radcor_R = np.array(radcorr_df[cut].R)

            for theta in binning["thetabins"]:
                # Cut data/mc for the w/q2/theta bin we're in
                data = rec[make_cuts(rec, w, q2, theta)].copy()
                data_e = empty[make_cuts(empty, w, q2, theta)].copy()
                data_mc = mc_rec[make_cuts(mc_rec, w, q2, theta)].copy()
                thrown = mc_thrown[make_cuts(mc_thrown, w, q2, theta)].copy()
                num_good = np.sum(data.cut_fid)
                if num_good < 6:
                    continue
                # elif num_good <= 24:
                #     bins = 24
                # else:
                #     bins = 10

                # Passes the plot section  at least once
                pass_plotting = True

                fig = plt.figure(
                    figsize=(12, 9), constrained_layout=True)
                gs = fig.add_gridspec(2, 1, height_ratios=[2, 1])
                ax1 = fig.add_subplot(gs[0])
                ax2 = fig.add_subplot(gs[1], sharex=ax1)
                maxs = None
                if overlap is not None:
                    for k, v in overlapSettings.items():
                        old_data = overlap_df[(overlap_df.W_min == w.left)
                                              & (overlap_df.Q2_min == q2.left)
                                              & (overlap_df.cos_t == theta.left)
                                              & (overlap_df.experiment == k)]

                        ebar = ax1.errorbar(old_data.phi, old_data.y, yerr=old_data.yerr,
                                            marker=v['symbol'], linestyle="",
                                            zorder=1, label=f"",
                                            markersize=10, alpha=0.4, c=v['color'])
                        maxs = np.max(old_data.y)*1.5
                        ct_ax[theta.left].errorbar(old_data.phi, old_data.y, yerr=old_data.yerr,
                                                   marker=v['symbol'], linestyle="",
                                                   markersize=5, alpha=0.4, c=v['color'])

                plot_maid_model(ax1, w, q2, theta, xs)
                maid_top = plot_maid_model(ct_ax[theta.left], w, q2, theta, xs)

                binCenter = binCetnerCorrection(w, q2, theta, num_bins=bins)

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
                        _empty_y, _ = hist_data(
                            data_e[data_e.cut_fid], density=False, bins=bins)
                        _mc_rec_y, _ = hist_data(
                            data_mc[data_mc.cut_fid], density=False, bins=bins)
                    elif cuts == 1:
                        # Histogram the data for plotting
                        marker = "^"
                        _data_y, _x = hist_data(
                            data[~data.cut_fid], density=False, bins=bins)
                        _empty_y, _ = hist_data(
                            data_e[~data_e.cut_fid], density=False, bins=bins)
                        _mc_rec_y, _ = hist_data(
                            data_mc[~data_mc.cut_fid], density=False, bins=bins)
                    else:
                        marker = 'd'
                        _data_y, _x = hist_data(
                            data, density=False, bins=bins)
                        _empty_y, _ = hist_data(
                            data_e, density=False, bins=bins)
                        _mc_rec_y, _ = hist_data(
                            data_mc, density=False, bins=bins)

                    _thrown_y, _ = hist_data(thrown, density=False, bins=bins)
                    # _data_y = (_data_y - _empty_y/Q_EMPTY)

                    # Remove points with 0 data count
                    cut = ~(_data_y == 0)
                    x = _x[cut]
                    mc_rec_y = _mc_rec_y[cut]
                    thrown_y = _thrown_y[cut]
                    data_y = _data_y[cut]
                    N_y = _data_y[cut]
                    N_empty = _empty_y[cut]

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

                    # Calculate acceptance and correct data
                    flux = virtual_photon_flux(w.left, q2.left) * luminosity()
                    stat_error = statistical(
                        N_y, N_empty, kin_bin_width, acceptance, flux)

                    # Normalize with bin widths
                    try:
                        data_y = data_y/kin_bin_width
                    except ValueError:
                        continue

                    y = data_y / acceptance / flux / binCenter(x) / radcor_R

                    error_bar = get_error_bars(
                        y, mc_rec_y, thrown_y, stat_error)

                    ebar = ax1.errorbar(x, y, yerr=error_bar,
                                        marker=marker, linestyle="",
                                        zorder=1, label=f"{name}", markersize=10, alpha=0.4)

                    out = fit_model(ax1, model_new, x, y, xs,
                                    ebar[0].get_color(), name)
                    ax1.legend(loc='upper right')
                    ax1.set_ylabel(
                        r'$\frac{\mathbf{d}\sigma}{\mathbf{d} \omega} \left[\frac{\mu b}{sr}\right]$')
                    ax1.set_xlabel(r'$\phi_{\pi}^{*}$')

                    ct_ax[theta.left].errorbar(x, y, yerr=error_bar,
                                               marker=marker, linestyle="",
                                               zorder=1,
                                               label=f"{plot_label[theta.left]}",
                                               markersize=5, alpha=0.8)

                    if cuts == 0:
                        ct_ax[theta.left].set_ylabel(
                            r'$\frac{\mathbf{d}\sigma}{\mathbf{d} \omega} \left[\frac{\mu b}{sr}\right]$')
                        ct_ax[theta.left].set_xlabel(r'$\phi_{\pi}^{*}$')
                        ct_ax[theta.left].legend(loc='upper right')
                        out = fit_model(ct_ax[theta.left], model_new, x, y, xs,
                                        ebar[0].get_color(), "")

                        try:
                            top = np.max(y)*1.5
                        except ValueError:
                            top = np.nan
                        if np.isnan(top) or np.isinf(top):
                            ax1.set_ylim(bottom=0, top=1.0)
                        else:
                            ax1.set_ylim(bottom=0, top=max(
                                top, max(maid_top, maxs)))
                            ct_ax[theta.left].set_ylim(
                                bottom=0.0, top=max(top, maid_top))

                ax1.set_title(f"$W$ : {w} , $Q^2$ : {q2} $\\theta$ : {theta}")
                fig.savefig(f"{out_folder}/crossSections/w_{w.left:0.3f}_q2_{q2.left:0.3f}_theta_{theta.left}.png",
                            bbox_inches='tight')

            if pass_plotting:
                CosTfig.suptitle(
                    f'$W~~[{w.left:0.3f},{w.right:0.3f})~~~~Q^2~~[{q2.left:0.3f}, {q2.right:0.3f})$', fontsize=16)
                CosTfig.align_ylabels()
                CosTfig.savefig(f"{out_folder}/crossSections/cost_w_{w.left:0.3f}_q2_{q2.left:0.3f}.png",
                                bbox_inches='tight', dpi=250)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Make Cross Sections")
    parser.add_argument("--mc", dest="mc_data_file_path",
                        type=str, help="MC csv file", required=True)
    parser.add_argument("--data", dest="rec_data_file_path",
                        type=str, help="Data csv file", required=True)
    parser.add_argument("--empty", dest="empty_file_path",
                        type=str, help="Empty run csv file", required=True)
    parser.add_argument("--folder", dest="out_folder", type=str,
                        help="Folder for plots", required=False, default="plots")
    parser.add_argument("--overlap", dest="overlap", type=str, help="Location of overlap data csv", required=False,
                        default=None)
    parser.add_argument("--radcorr", dest="radcorr", type=str, help="Location of radcorr data csv", required=False,
                        default=None)
    args = parser.parse_args()

    # Start to main

    print("Start setup")
    start = time.time_ns()
    # Load mc file
    mc_rec, mc_thrown = read_csv(args.mc_data_file_path)
    # Load reconstructed file
    rec = read_csv(args.rec_data_file_path, True)
    empty_target = read_csv(args.empty_file_path, True)

    # Cut for missing mass
    rec, mc_rec, empty_target = cut_for_MM(rec, mc_rec, empty_target)

    # Make bins in the dataframes from the bins above
    rec = prep_for_ana(rec, w_bins, q2_bins, theta_bins)
    empty_target = prep_for_ana(empty_target, w_bins, q2_bins, theta_bins)
    mc_rec = prep_for_ana(mc_rec, w_bins, q2_bins, theta_bins)
    mc_thrown = prep_for_ana(mc_thrown, w_bins, q2_bins, theta_bins)

    # Create dict of sorted bins
    binning = dict()
    binning["wbins"] = pd.Index.sort_values(pd.unique(rec.w_bin))
    binning["q2bins"] = pd.Index.sort_values(pd.unique(rec.q2_bin))
    binning["thetabins"] = pd.Index.sort_values(pd.unique(rec.theta_bin))
    end = time.time_ns()
    print(f"Done setup: {(end-start)/1E9:0.2f}Sec")

    main(rec, mc_rec, mc_thrown, empty_target, binning, bins=12,
         out_folder=args.out_folder, overlap=args.overlap, radcorr=args.radcorr)
