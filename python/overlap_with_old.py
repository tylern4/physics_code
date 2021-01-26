#!/usr/bin/env python

from typing import Dict
import matplotlib  # noqa

matplotlib.use("agg")  # noqa
import warnings  # noqa

warnings.filterwarnings("ignore")  # noqa

import datetime
import argparse
import inspect
import os
import time
import matplotlib.pyplot as plt
from tqdm import tqdm
from calc_xsections import *


def draw_cos_bin(data, mc_rec_data, thrown_data, w, q2, cos_t_bins, out_folder, bins, models_fits={"model": model_new}):
    df = pd.read_csv(
        "/Users/tylern/github/crossSectionPlotter/overlap_data.csv", index_col=0)

    which_plot = {
        -1.0: [0, 0],
        -0.8: [0, 1],
        -0.6: [1, 0],
        -0.4: [1, 1],
        -0.2: [2, 0],
        0.0: [2, 1],
        0.2: [3, 0],
        0.4: [3, 1],
        0.6: [4, 0],
        0.8: [4, 1]
    }
    plot_label = {
        -1.0: "$\cos(\\theta)=(-1.0,-0.8]$",
        -0.8: "$\cos(\\theta)=(-0.8,-0.6]$",
        -0.6: "$\cos(\\theta)=(-0.6,-0.4]$",
        -0.4: "$\cos(\\theta)=(-0.4,-0.2]$",
        -0.2: "$\cos(\\theta)=(-0.2,0.0]$",
        0.0: "$\cos(\\theta)=(0.0,0.2]$",
        0.2: "$\cos(\\theta)=(0.2,0.4]$",
        0.4: "$\cos(\\theta)=(0.4,0.6]$",
        0.6: "$\cos(\\theta)=(0.6,0.8]$",
        0.8: "$\cos(\\theta)=(0.8,1.0]$"
    }

    fig, ax = plt.subplots(5, 2, figsize=(
        12, 9), sharex=True, gridspec_kw={'hspace': 0.1})

    fig.suptitle(
        f"W=({np.round(w.left,3)}, {np.round(w.right,3)}],\t$Q^2$=({np.round(q2.left,3)}, {np.round(q2.right,3)}]", fontsize=20
    )
    color = '0.5'
    alpha = 0.3
    sizes = fig.get_size_inches()
    size = 8.0 * np.sqrt(sizes[0]**2.0 + sizes[1]**2.0)
    fig.text(0.5, 0.5, 'Preliminary',
             fontsize=size, color=color,
             ha='center', va='center', alpha=alpha, rotation=45, zorder=0)

    phi_bins = np.linspace(0, 2 * np.pi, 200)
    thetabins = pd.unique(cos_t_bins)

    for c, cos_t in enumerate(thetabins):
        old_data = df[(df.W_min == w.left) & (
            df.Q2_min == q2.left) & (df.cos_t == cos_t.left)]

        a = which_plot[round(cos_t.left, 1)][0]
        b = which_plot[round(cos_t.left, 1)][1]
        cos_label = plot_label[round(cos_t.left, 1)]

        _w = (w.left + w.right) / 2.0
        _q2 = (q2.left + q2.right) / 2.0
        _cos_t = (cos_t.left + cos_t.right) / 2.0

        _data = data[cos_t == data.theta_bin].copy()
        _mc_rec_data = mc_rec_data[cos_t == mc_rec_data.theta_bin].copy()
        _thrown_data = thrown_data[cos_t == thrown_data.theta_bin].copy()

        # flux = virtual_photon(w.left, q2.left, ENERGY)
        flux = virtual_photon_flux(w.left, q2.left) * luminosity()

        xs = np.linspace(0, 2 * np.pi, 100)
        # print(len(_data.cut_fid.to_numpy()), len(_data.phi.to_numpy()))
        cut_fids = {"With Fid cuts": True,
                    "No fid cuts": False
                    }
        for name, cuts in cut_fids.items():
            if cuts:
                # Histogram the data for plotting
                data_y, x = hist_data(
                    _data[_data.cut_fid], density=False, bins=bins)
                mc_rec_y, _ = hist_data(
                    _mc_rec_data[_mc_rec_data.cut_fid], density=False, bins=bins)
                _ax = ax[a][b].twinx()
                _ax.errorbar(
                    old_data.phi,
                    old_data.y,
                    yerr=old_data.yerr,
                    marker=".",
                    linestyle="",
                    c='r',
                    zorder=1,
                    label=f"E-99-107" if c == 0 else None,
                )
                plot_maid_model(_ax, w, q2, cos_t, xs)
            else:
                # Histogram the data for plotting
                data_y, x = hist_data(data, density=False, bins=bins)
                mc_rec_y, _ = hist_data(_mc_rec_data, density=False, bins=bins)
                # _ax.plot(phis, crossSections, c='r', linestyle='dotted')
                # _ax.set_ylim(bottom=0, top=np.max(crossSections)*1.5)
                # _ax.text(0.0, 1.2*np.max(crossSections), cos_label)

            thrown_y, _ = hist_data(_thrown_data, density=True, bins=bins)

            # Remove points with 0 data count
            cut = ~(data_y == 0)
            x = x[cut]
            mc_rec_y = mc_rec_y[cut]
            thrown_y = thrown_y[cut]
            data_y = data_y[cut]

            # Remove points with 0 mc_rec count and put 1???
            mc_rec_y = np.where(mc_rec_y == 0, 1, mc_rec_y)
            thrown_y = np.where(thrown_y == 0, 1, thrown_y)
            acceptance = mc_rec_y / thrown_y

            delta_W = (w.right-w.left)
            delta_Q2 = (q2.right-q2.left)
            delta_Theta = 0.2
            __phis = np.linspace(0, 2 * np.pi, bins)
            delta_phi = __phis[1] - __phis[0]
            kin_bin_width = delta_W * delta_Q2 * delta_Theta * delta_phi

            try:
                data_y = data_y / kin_bin_width
            except ValueError as e:
                print(e)
                continue

            flux = virtual_photon_flux(w.left, q2.left) * luminosity()
            y = data_y / acceptance / flux
            error_bar = get_error_bars(y, mc_rec_y, thrown_y)

            # error_bar = np.ones_like(y) * 0.1

            F = (mc_rec_y/thrown_y)
            error = np.sqrt(((thrown_y-mc_rec_y)*mc_rec_y) /
                            np.power(thrown_y, 3))/F
            error_bar = np.sqrt(
                np.power((y*error), 2) + np.power(stats.sem(y), 2))

            error_bar = stats.sem(acceptance)

            try:
                ax[a][b].set_ylim(bottom=0, top=np.max(y)*1.5)
            except ValueError as e:
                print(e)

            ebar = ax[a][b].errorbar(
                x,
                y,
                yerr=error_bar,
                marker=".",
                linestyle="",
                # c="k",
                zorder=1,
                label=f"{name}" if c == 0 else None,
            )

            # Try dropping 0's before fitting
            if y.size <= 5:
                continue
            for name, func in models_fits.items():
                # Make model from function
                model = Model(func)
                # Make fit parameters
                params = model.make_params()
                # Make sure to set inital values to 1
                # so fit doesn't fail
                for p in params:
                    params[p].set(value=1)

                # Fit the model
                try:
                    out = model.fit(y, params, x=x)
                    ax[a][b].plot(xs, out.eval(params=out.params, x=xs),
                                  linewidth=2.0, c=ebar[0].get_color(), alpha=0.1)
                except ValueError as e:
                    print(e)
                # Plot the fitted model with output parameters and same x's as model

                # # Get uncertinty and plot between 2 sigmas
                # dely = out.eval_uncertainty(sigma=2, x=xs)
                # ax[a][b].fill_between(xs, out.eval(params=out.params, x=xs)-dely,
                #                       out.eval(
                #     params=out.params, x=xs)+dely,
                #     color=ebar[0].get_color(), alpha=0.1,
                #     label='2-$\sigma$ uncertainty band')

    if not os.path.exists(f'{out_folder}/CosT'):
        os.makedirs(f'{out_folder}/CosT')
    fig.legend(loc="upper right")
    fig.savefig(
        f"{out_folder}/CosT/W[{np.round(w.left,3)},{np.round(w.right,3)}]_Q2[{np.round(q2.left,3)},{np.round(q2.right,3)}]_{bins}_CosT.png", bbox_inches='tight'
    )


def draw_cos_plots(data, mc_rec_data, thrown_data, w, q2, cos_t_bins, out_folder):
    # with multiprocessing.Pool() as pool:
    #     inputs = []
    #     for bins in range(10, 11):
    #         inputs.append((func, data, mc_rec_data, thrown_data,
    #                        w, q2, cos_t_bins, out_folder, bins))
    #     pool.starmap(draw_cos_bin, inputs)
    # for bins in range(10, 26, 2):
    for bins in [24]:
        draw_cos_bin(data, mc_rec_data, thrown_data,
                     w, q2, cos_t_bins, out_folder, bins, models_fits={"new": model_new})

    # draw_cos_bin(data, mc_rec_data, thrown_data,
    #              w, q2, cos_t_bins, out_folder, 10, models_fits={"new": model_new})


def draw_xsec_plots(func, data, mc_rec_data, thrown_data, w, q2, cos_t, out_folder, bins):
    fig, ax = plt.subplots(2, 2, figsize=(12, 9))
    fig.suptitle(
        f"W={w},\t$Q^2$={q2},\tcos($\\theta$)={cos_t} \n{bins} bins in $\phi$"
    )
    xs = np.linspace(0, 2 * np.pi, 100)
    data_y, data_x = bh.numpy.histogram(
        data.phi, bins=bins, range=(0, 2 * np.pi), threads=4)
    x = (data_x[1:] + data_x[:-1]) / 2.0
    mc_rec_y, _ = bh.numpy.histogram(
        mc_rec_data.phi, bins=bins, range=(0, 2 * np.pi), threads=4
    )
    thrown_y, _ = bh.numpy.histogram(
        thrown_data.phi, bins=bins, range=(0, 2 * np.pi), threads=4
    )

    # Change 0's to 1 for division
    thrown_y = np.where(thrown_y == 0, 1, thrown_y)
    mc_rec_y = np.where(mc_rec_y == 0, 1, mc_rec_y)

    # stats.sem(thrown_y)
    ax[0][0].errorbar(
        x,
        thrown_y,
        marker=".",
        yerr=stats.sem(thrown_y),
        c="r",
        linestyle="",
        label="thrown",
    )
    _ax = ax[0][0].twinx()
    _ax.errorbar(
        x,
        mc_rec_y,
        marker=".",
        yerr=stats.sem(mc_rec_y),
        c="orange",
        linestyle="",
        label="mc_rec",
    )
    ax[0][1].errorbar(
        x, data_y, yerr=stats.sem(data_y), marker=".", linestyle="", label="data",
    )

    acceptance = np.nan_to_num(thrown_y / mc_rec_y)
    ax[1][0].errorbar(
        x,
        1/acceptance,
        yerr=0,
        marker=".",
        c="g",
        linestyle="",
        label="acceptance",
    )

    y = data_y * acceptance
    y = y / np.max(y)

    error_bar = np.ones_like(y)

    F = mc_rec_y/thrown_y
    error = np.sqrt(((thrown_y-mc_rec_y)*mc_rec_y) /
                    np.power(thrown_y, 3))/F
    error_bar = np.sqrt(
        np.power((y*error), 2) + np.power(stats.sem(y), 2))

    ax[1][1].errorbar(
        x,
        y,
        yerr=error_bar,
        marker=".",
        linestyle="",
        c="k",
        zorder=1,
        label="corrected",
    )
    ax[1][1].set_ylim(bottom=0, top=np.max(y)*1.5)

    phi_bins = np.linspace(0, 2 * np.pi, 200)
    crossSections = []
    phis = []

    _w = (w.left + w.right) / 2.0
    _q2 = (q2.left + q2.right) / 2.0
    _cos_t = (cos_t.left + cos_t.right) / 2.0

    for phi in phi_bins:
        crossSections.append(
            maid(ENERGY, _w, _q2, _cos_t, np.degrees(phi)))
        phis.append(phi)

    crossSections = np.array(crossSections)
    phis = np.array(phis)
    _ax = ax[1][1].twinx()
    _ax.plot(phis, crossSections, c='r',
             label='maid2007', linestyle='dotted')

    _ax.set_ylim(bottom=0, top=np.max(crossSections)*1.5)

    # popt, pcov = curve_fit(func, x, y, maxfev=8000)
    # To compute one standard deviation errors on the parameters use
    # https://stackoverflow.com/questions/49130343/is-there-a-way-to-get-the-error-in-fitting-parameters-from-scipy-stats-norm-fit
    # perr = np.sqrt(np.diag(pcov))
    # ax[1][1].plot(xs, func(xs, *popt), c="#9467bd",
    #               linewidth=2.0)

    fig.legend()
    if not os.path.exists(f'{out_folder}/acc'):
        os.makedirs(f'{out_folder}/acc')

    plt.savefig(
        f"{out_folder}/acc/W[{w.left},{w.right}]_Q2[{q2.left},{q2.right}]_cos(theta)[{cos_t.left},{cos_t.right}]_{bins}.png"
    )


def draw_plots(func, data, mc_rec_data, thrown_data, w, q2, cos_t, out_folder):
    # with multiprocessing.Pool() as pool:
    #     inputs = []
    #     for bins in range(10, 11):
    #         inputs.append((func, data, mc_rec_data, thrown_data,
    #                        w, q2, cos_t, out_folder, bins))
    #     pool.starmap(draw_xsec_plots, inputs)

    draw_xsec_plots(func, data, mc_rec_data, thrown_data,
                    w, q2, cos_t, out_folder, 10)


def draw_xsection(rec: pd.DataFrame, mc_rec: pd.DataFrame, thrown: pd.DataFrame, func, out_folder: str, binning: Dict):
    # executor = get_reusable_executor(max_workers=len(np.unique(rec.theta_bin)))
    total_num = (
        len(binning["wbins"])
        * len(binning["q2bins"])
    )

    pbar = tqdm(total=total_num)

    for w in binning["wbins"]:
        for q2 in binning["q2bins"]:
            #################################################
            rec_cut = (
                (w == rec.w_bin) & (q2 == rec.q2_bin)
            )
            mc_rec_cut = (
                (w == mc_rec.w_bin)
                & (q2 == mc_rec.q2_bin)
            )
            thrown_cut = (
                (w == thrown.w_bin)
                & (q2 == thrown.q2_bin)
            )

            data = rec[rec_cut]
            mc_rec_data = mc_rec[mc_rec_cut]
            thrown_data = thrown[thrown_cut]
            draw_cos_plots(data, mc_rec_data,
                           thrown_data, w, q2, rec.theta_bin, out_folder)
            #################################################
            # for cos_t in thetabins:
            #     rec_cut = (
            #         (w == rec.w_bin) & (q2 == rec.q2_bin) & (
            #             cos_t == rec.theta_bin)
            #     )
            #     mc_rec_cut = (
            #         (w == mc_rec.w_bin)
            #         & (q2 == mc_rec.q2_bin)
            #         & (cos_t == mc_rec.theta_bin)
            #     )
            #     thrown_cut = (
            #         (w == thrown.w_bin)
            #         & (q2 == thrown.q2_bin)
            #         & (cos_t == thrown.theta_bin)
            #     )

            #     data = rec[rec_cut]
            #     mc_rec_data = mc_rec[mc_rec_cut]
            #     thrown_data = thrown[thrown_cut]

            #     # draw_plots(
            #     #     func, data, mc_rec_data, thrown_data, w, q2, cos_t, out_folder
            #     # )
            pbar.update(1)

    pbar.close()


def draw_kinematics(rec, w_bins, q2_bins, theta_bins, name="reconstructed"):
    rec = rec.dropna()

    fig, ax = plt.subplots(figsize=(12, 9))
    h = ax.hist2d(rec.w.to_numpy(), rec.q2.to_numpy(),
                  bins=200,
                  range=[[np.min(w_bins), np.max(w_bins)],
                         [np.min(q2_bins), np.max(q2_bins)]],
                  # cmin=1
                  )

    for w in w_bins:
        ax.axvline(w, c='w')
    for q2 in q2_bins:
        ax.axhline(q2, c='w')
    ax.set_xlabel(f"W [$\mathrm{{{{GeV}}}}$]", fontsize=30)
    ax.set_ylabel(f"$Q^2$ [$\mathrm{{{{GeV}}}}^2$]", fontsize=30)
    ax.tick_params(axis='both', which='major', labelsize=15)
    fig.suptitle(f"W vs $Q^2$", fontsize=30)
    fig.colorbar(h[3], ax=ax)

    if not os.path.exists(f'{out_folder}/kinematics'):
        os.makedirs(f'{out_folder}/kinematics')

    fig.savefig(f"{out_folder}/kinematics/W_vs_Q2_{name}.png",
                bbox_inches='tight')

    #########################################
    fig2, ax2 = plt.subplots(figsize=(12, 9))
    y, x = bh.numpy.histogram(rec.w.to_numpy(), bins=250, range=[
        np.min(w_bins), np.max(w_bins)])
    x = (x[1:] + x[:-1])/2.0

    ax2.errorbar(x, y, yerr=stats.sem(y), linestyle="", marker=".")
    ax2.set_xlabel(f"W [$GeV/c^2$]", fontsize=30)
    ax2.tick_params(axis='both', which='major', labelsize=15)
    ax2.set_xlim(np.min(w_bins), np.max(w_bins))
    ax2.set_ylim(bottom=0)

    fig2.suptitle(f"W $(\mathrm{{{{n}}}}\pi^+)$", fontsize=30)
    fig2.savefig(f"{out_folder}/kinematics/W_{name}.png", bbox_inches='tight')
    #########################################

    fig1, ax1 = plt.subplots(figsize=(12, 9))
    h = ax1.hist2d(rec.cos_theta.to_numpy(),
                   rec.phi.to_numpy(), bins=100)
    for t in theta_bins:
        ax1.axvline(t, c='w')

    for phi in np.linspace(0, 2*np.pi, 10):
        ax1.axhline(phi, c='w')

    ax1.set_xlabel(f"$\cos(\\theta)$", fontsize=30)
    ax1.set_ylabel(f"$\\phi$", fontsize=30)
    ax1.tick_params(axis='both', which='major', labelsize=15)
    fig1.suptitle(f"$\cos(\\theta)$ vs $\\phi$", fontsize=30)
    fig1.colorbar(h[3], ax=ax1)

    fig1.savefig(f"{out_folder}/kinematics/theta_vs_phi_{name}.png",
                 bbox_inches='tight')
    #########################################
    fig3, ax3 = plt.subplots(ncols=1, nrows=2, figsize=(
        16, 16), sharex=True, gridspec_kw={'hspace': 0.05})
    h = ax3[0].hist2d(rec.w.to_numpy(), rec.q2.to_numpy(),
                      bins=200, range=[[np.min(w_bins), np.max(w_bins)], [np.min(q2_bins), np.max(q2_bins)]])

    ax3[1].errorbar(x, y, yerr=stats.sem(y), linestyle="", marker=".", ms=10)
    ax3[1].set_xlabel(f"W [$\mathrm{{{{GeV}}}}$]", fontsize=30)
    ax3[0].set_ylabel(f"$Q^2$ [$\mathrm{{{{GeV}}}}^2$]", fontsize=30)
    ax3[0].tick_params(axis='y', which='major', labelsize=25)
    ax3[0].tick_params(axis='x', which='major', labelsize=0)
    ax3[1].tick_params(axis='both', which='major', labelsize=25)
    ax3[1].set_xlim(np.min(w_bins), np.max(w_bins))
    for axxx in ax3:
        axxx.label_outer()
    ax3[0].set_title(f"W vs $Q^2$", fontsize=30)
    fig3.savefig(f"{out_folder}/kinematics/wq2_{name}.png",
                 bbox_inches='tight')

    fig4, ax4 = plt.subplots(figsize=(12, 9))
    H, xedges, yedges = bh.numpy.histogram2d(
        rec.w.to_numpy(), rec.q2.to_numpy(), bins=(w_bins, q2_bins))
    H = H.T  # Let each row list bins with common y range.
    X, Y = np.meshgrid(xedges, yedges)
    im = ax4.pcolormesh(X, Y, H)
    fig4.colorbar(im, ax=ax4)
    ax4.set_title(f"W vs $Q^2$", fontsize=30)
    ax4.set_xlabel(f"W [$\mathrm{{{{GeV}}}}$]", fontsize=30)
    ax4.set_ylabel(f"$Q^2$ [$\mathrm{{{{GeV}}}}^2$]", fontsize=30)
    fig4.savefig(f"{out_folder}/kinematics/wq2_binned_{name}.png",
                 bbox_inches='tight')

    fig4, ax4 = plt.subplots(figsize=(12, 9))
    try:
        H, xedges, yedges = bh.numpy.histogram2d(
            rec.cos_theta.to_numpy(),
            rec.phi.to_numpy(),
            bins=(11, 10))
    except ValueError as ve:
        print(rec.cos_theta.to_numpy().size)
        print(rec.phi.to_numpy().size)
        return

    H = H.T  # Let each row list bins with common y range.
    X, Y = np.meshgrid(xedges, yedges)
    im = ax4.pcolormesh(X, Y, H)
    fig4.colorbar(im, ax=ax4)
    ax4.set_xlabel(f"$\cos(\\theta)$", fontsize=30)
    ax4.set_ylabel(f"$\\phi$", fontsize=30)
    ax4.tick_params(axis='both', which='major', labelsize=15)
    fig1.suptitle(f"$\cos(\\theta)$ vs $\\phi$", fontsize=30)
    fig4.savefig(f"{out_folder}/kinematics/cos_phi_binned_{name}.png",
                 bbox_inches='tight')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Make Cross Sections")
    parser.add_argument(
        "--mc", dest="mc_data_file_path", type=str, help="MC csv file", required=True
    )
    parser.add_argument(
        "--data",
        dest="rec_data_file_path",
        type=str,
        help="Data csv file",
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
    parser.add_argument(
        "--draw_kin",
        action='store_true'
    )
    parser.add_argument(
        "--e1f",
        action='store_true'
    )

    args = parser.parse_args()
    if args.e1f:
        ENERGY = 5.479

    total_time = time.time()
    mc_data_file_path = args.mc_data_file_path
    rec_data_file_path = args.rec_data_file_path
    out_folder = args.out_folder

    start = time.time()
    mc_rec, mc_thrown = read_csv(mc_data_file_path)
    stop = time.time()
    # print(f"\n\nread time mc_df: {stop - start}\n\n")
    # print(f"\n\ntime: {stop - total_time}\n\n")

    # mc_rec = mc_rec[(mc_rec.w > 0) & (mc_rec.mm2 > 0.5) & (mc_rec.mm2 < 1.5)]
    mc_rec["cos_theta"] = np.cos(mc_rec.theta).astype(np.float32)

    mc_thrown = mc_thrown[(mc_thrown.w > 0)]
    mc_thrown["cos_theta"] = np.cos(mc_thrown.theta).astype(np.float32)

    start = time.time()
    # rec = feather.read_feather(rec_data_file_path)
    rec = read_csv(rec_data_file_path, True)
    stop = time.time()
    # print(f"\n\nread time rec: {stop - start}\n\n")
    # rec = rec[(rec.w > 0) & (rec.mm2 > 0.5) & (rec.mm2 < 1.5)]
    rec["cos_theta"] = np.cos(rec.theta).astype(np.float32)

    sector_cuts = mm_cut(rec, sigma=10)

    cuts = False
    mc_cuts = False

    for sec, min_max in sector_cuts.items():
        cuts |= (
            (rec.electron_sector == sec)
            & (rec.mm2 >= min_max[0])
            & (rec.mm2 <= min_max[1])
        )
        mc_cuts |= (
            (mc_rec.electron_sector == sec)
            & (mc_rec.mm2 >= min_max[0])
            & (mc_rec.mm2 <= min_max[1])
        )
    rec = rec[cuts]
    mc_rec = mc_rec[mc_cuts]
    mc_rec = mc_rec[["w", "q2", "mm2", "cos_theta",
                     "phi", "helicty", "electron_sector", "cut_fid"]].copy(deep=True)
    mc_thrown = mc_thrown[["w", "q2", "mm2", "cos_theta", "phi", "helicty", "electron_sector", "cut_fid"]].copy(
        deep=True
    )
    rec = rec[["w", "q2", "mm2", "cos_theta",
               "phi", "cut_fid", "helicty",  "electron_sector"]].copy(deep=True)

    # Specifically put in bin edges
    # TODO ##################### BINS ######################
    w_bins = np.array([1.1, 1.12, 1.14, 1.16, 1.18, 1.2, 1.22, 1.24, 1.26, 1.28, 1.3,
                       1.32, 1.34, 1.36, 1.38, 1.4, 1.42, 1.44, 1.46, 1.48, 1.5, 1.52,
                       1.54, 1.56, 1.58, 1.6, 1.62, 1.64, 1.66, 1.68, 1.7, 1.72, 1.74,
                       1.76, 1.78, 1.8])
    q2_bins = np.array([1.1, 1.30, 1.56, 1.87, 2.23, 2.66])
    theta_bins = np.array([-1.0, -0.8, -0.6, -0.4, -0.2,
                           0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2])
    # TODO ##################### BINS ######################
    if args.draw_kin:
        draw_kinematics(rec, w_bins, q2_bins, theta_bins)
        draw_kinematics(mc_rec, w_bins, q2_bins, theta_bins, "mc_rec")
        draw_kinematics(mc_thrown, w_bins, q2_bins, theta_bins, "thrown")

        for sec in range(1, 7):
            sec_data = rec[rec.electron_sector == sec]
            draw_kinematics(sec_data, w_bins, q2_bins,
                            theta_bins, f"rec_{sec}")

            sec_mc_rec = mc_rec[mc_rec.electron_sector == sec]
            draw_kinematics(sec_mc_rec, w_bins, q2_bins,
                            theta_bins, f"mc_rec_{sec}")

            sec_mc_thrown = mc_thrown[mc_thrown.electron_sector == sec]
            draw_kinematics(sec_mc_thrown, w_bins, q2_bins,
                            theta_bins, f"thrown_{sec}")

    mc_rec["w_bin"] = pd.cut(mc_rec["w"], bins=w_bins, include_lowest=False)
    mc_rec["q2_bin"] = pd.cut(mc_rec["q2"], bins=q2_bins, include_lowest=False)
    mc_rec["theta_bin"] = pd.cut(
        mc_rec["cos_theta"], bins=theta_bins, include_lowest=False
    )

    mc_thrown["w_bin"] = pd.cut(
        mc_thrown["w"], bins=w_bins, include_lowest=False)
    mc_thrown["q2_bin"] = pd.cut(
        mc_thrown["q2"], bins=q2_bins, include_lowest=False)
    mc_thrown["theta_bin"] = pd.cut(
        mc_thrown["cos_theta"], bins=theta_bins, include_lowest=False
    )

    rec["w_bin"] = pd.cut(rec["w"], bins=w_bins, include_lowest=False)
    rec["q2_bin"] = pd.cut(rec["q2"], bins=q2_bins, include_lowest=False)
    rec["theta_bin"] = pd.cut(
        rec["cos_theta"], bins=theta_bins, include_lowest=False)

    mc_rec.dropna(inplace=True)
    mc_thrown.dropna(inplace=True)
    rec.dropna(inplace=True)

    print(f"===========================\nmc_rec:\n\n")
    print(f"{mc_rec.info(verbose=True, memory_usage='deep')}")
    print(f"\n\n===========================")
    print(f"===========================\nmc_thrown:\n\n")
    print(f"{mc_thrown.info(verbose=True, memory_usage='deep')}")
    print(f"\n\n===========================")
    print(f"===========================\nrec:\n\n")
    print(f"{rec.info(verbose=True, memory_usage='deep')}")
    print(f"\n\n===========================")

    binning = dict()
    binning["wbins"] = pd.Index.sort_values(pd.unique(rec.w_bin))
    binning["q2bins"] = pd.Index.sort_values(pd.unique(rec.q2_bin))
    binning["thetabins"] = pd.Index.sort_values(pd.unique(rec.theta_bin))

    draw_xsection(rec, mc_rec, mc_thrown, model_new,
                  out_folder, binning)

    stop = time.time()
    print(f"\n\nFull Running time: {stop - total_time}\n\n")
    print(
        f"\n\nFull Running time: {datetime.timedelta(seconds=(time.time() - total_time))}\n\n"
    )
