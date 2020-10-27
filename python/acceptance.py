#!/usr/bin/env python

import matplotlib  # noqa
matplotlib.use("agg")  # noqa

import warnings
# from loky import get_reusable_executor
import multiprocessing
import os
from maid_interface import maid_2007_Npi as maid
import datetime
import boost_histogram as bh
import pyarrow
from pyarrow import feather
from pyarrow import csv
import pyarrow as pa
import time
import argparse
from scipy.optimize import curve_fit
from scipy import stats
from scipy.special import erfc
from tqdm import tqdm
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import numba
#from fit_functions import gauss, degauss, half_max_x, model

warnings.filterwarnings("ignore")

ENERGY = 4.81726


def read_csv(file_name):
    names = [
        "electron_sector",
        "w",
        "q2",
        "theta",
        "phi",
        "mm2",
        "helicty",
        "photon_flux",
        "type"
    ]
    dtype = {
        "electron_sector": "int8",
        "helicty": "int8",
        "w": "float32",
        "q2": "float32",
        "theta": "float32",
        "phi": "float32",
        "mm2": "float32",
        "photon_flux": "float32"
    }

    start = time.time()
    pyTable = csv.read_csv(
        file_name,
        read_options=csv.ReadOptions(use_threads=True, column_names=names),
        convert_options=csv.ConvertOptions(column_types=dtype),
    )
    # pyTable = pyTable.drop(["hash"])
    df = pyTable.to_pandas(strings_to_categorical=True)

    mc_rec = df[df.type == "mc_rec"]
    thrown = df[df.type == "thrown"]
    del df

    stop = time.time()
    # print(f"read_csv: {stop - start}")

    return (
        mc_rec,
        thrown,
    )


def model(x, a, b, c):
    """
    a => sigma_l + sigma_t
    b => epsilon*sigma_tt
    c => Sqrt(2epsilon(1+epsilon))* sigma_lt
    """
    f = a + b * np.cos(2 * x) + c * np.cos(x)
    return f


def degauss(x, A, mu, sigma, lambda1, lambda2):
    mu1 = sigma * sigma * lambda1 + x - mu
    mu2 = -sigma * sigma * lambda2 + x - mu
    ret = (
        A
        * 0.5
        / (1.0 / lambda1 + 1.0 / lambda2)
        * (
            np.exp(0.5 * np.power(sigma * lambda1, 2) + lambda1 * (x - mu))
            * erfc(mu1 / (sigma * np.sqrt(2.0)))
            + np.exp(0.5 * np.power(sigma * lambda2, 2) - lambda2 * (x - mu))
            * erfc(-mu2 / (sigma * np.sqrt(2.0)))
        )
    )
    return ret


def gauss(x, A, mu, sig):
    ret = np.exp(-np.power(x - mu, 2.0) / (2 * np.power(sig, 2.0)))
    return A * ret


def peak(x, c):
    return np.exp(-np.power(x - c, 2) / 16.0)


def lin_interp(x, y, i, half):
    return x[i] + (x[i + 1] - x[i]) * ((half - y[i]) / (y[i + 1] - y[i]))


def half_max_x(x, y):
    half = np.max(y) / 2.0
    signs = np.sign(np.add(y, -half))
    zero_crossings = signs[0:-2] != signs[1:-1]
    zero_crossings_i = np.where(zero_crossings)[0]
    return [
        lin_interp(x, y, zero_crossings_i[0], half),
        lin_interp(x, y, zero_crossings_i[1], half),
    ]


def mm_cut(df, sigma=6):
    data = {}
    fig, ax = plt.subplots(2, 3, figsize=(
        12, 9), sharex=True, sharey=True, gridspec_kw={'hspace': 0.0, 'wspace': 0.0})
    which_plot = {
        1: [0, 0],
        2: [0, 1],
        3: [0, 2],
        4: [1, 0],
        5: [1, 1],
        6: [1, 2],
    }
    for sec in range(1, 7):
        a = which_plot[sec][0]
        b = which_plot[sec][1]
        plt.figure(figsize=(12, 9))
        y, x = bh.numpy.histogram(
            df[df.electron_sector == sec].mm2, bins=500, density=True
        )
        x = (x[1:] + x[:-1]) / 2
        popt_g, pcov_g = curve_fit(gauss, x, y, maxfev=8000)
        #plt.plot(x, gauss(x, *popt_g), linewidth=2.0)
        plt.errorbar(x, y, yerr=stats.sem(y), fmt=".", zorder=1)
        ax[a][b].errorbar(
            x, y, yerr=stats.sem(y), fmt=".", zorder=1)

        # plt.axvline(popt_g[1] + sigma * popt_g[2])
        # plt.axvline(popt_g[1] - sigma * popt_g[2])

        p0 = [popt_g[0], popt_g[1], popt_g[2], 1.0, 1.0]
        popt, pcov = curve_fit(degauss, x, y, maxfev=8000)

        plt.plot(x, degauss(x, *popt), c="#9467bd", linewidth=2.0)
        ax[a][b].plot(x, degauss(
            x, *popt), c="#9467bd", linewidth=3.0, label=f"Sector {sec}")
        ax[a][b].legend()
        ax[a][b].set_xlabel(f"Mass $[ $\mathrm{{{{GeV}}}}^2 ]$")

        # find the FWHM
        xs = np.linspace(0.7, 1.5, 100000)
        hmx = half_max_x(xs, degauss(xs, *popt))
        fwhm = hmx[1] - hmx[0]
        plt.axvline(popt[1] + sigma * fwhm / 2.355, c="#9467bd")
        plt.axvline(popt[1] - sigma * fwhm / 2.355, c="#9467bd")

        ax[a][b].axvline(
            popt[1] + sigma * fwhm / 2.355, c="#9467bd", linewidth=3.0)
        ax[a][b].axvline(
            popt[1] - sigma * fwhm / 2.355, c="#9467bd", linewidth=3.0)

        if not os.path.exists(f'{out_folder}/cuts'):
            os.makedirs(f'{out_folder}/cuts')

        plt.title(
            f"Missing Mass Squared $e\left( p, \pi^{{{'+'}}} X \\right)$ in sector {sec}")

        plt.savefig(f"{out_folder}/cuts/MM2_cut_{sec}.png",
                    bbox_inches='tight')

        data[sec] = (popt[1] - sigma * fwhm / 2.355,
                     popt[1] + sigma * fwhm / 2.355)
    fig.suptitle(
        f"Missing Mass Squared $e\left( p, \pi^{{{'+'}}} X \\right)$", fontsize=20)
    fig.savefig(f"{out_folder}/cuts/MM2_cut_all.png",
                bbox_inches='tight')
    return data


def draw_cos_bin(func, data, mc_rec_data, thrown_data, w, q2, cos_t_bins, out_folder, bins):
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
    if str(np.round(q2.left, 3)) == "0.999":
        q2_left = 1.0
    else:
        q2_left = np.round(q2.left, 3)

    fig, ax = plt.subplots(5, 2, figsize=(
        12, 9), sharex=True, gridspec_kw={'hspace': 0.1})

    fig.suptitle(
        f"W=({np.round(w.left,3)}, {np.round(w.right,3)}],\t$Q^2$=({q2_left}, {np.round(q2.right,3)}]", fontsize=20
    )
    color = '0.5'
    alpha = 0.3
    sizes = fig.get_size_inches()
    size = 8.0 * np.sqrt(sizes[0]**2.0 + sizes[1]**2.0)
    fig.text(0.5, 0.5, 'Preliminary',
             fontsize=size, color=color,
             ha='center', va='center', alpha=alpha, rotation=45, zorder=0)

    thetabins = pd.unique(cos_t_bins)
    q2_left = 0.0
    for c, cos_t in enumerate(thetabins):
        a = which_plot[round(cos_t.left, 1)][0]
        b = which_plot[round(cos_t.left, 1)][1]
        cos_label = plot_label[round(cos_t.left, 1)]

        _data = data[cos_t == data.theta_bin]
        _mc_rec_data = mc_rec_data[cos_t == mc_rec_data.theta_bin]
        _thrown_data = thrown_data[cos_t == thrown_data.theta_bin]

        flux = np.sum(_data['photon_flux'].to_numpy()) / \
            len(_data['photon_flux'].to_numpy())

        xs = np.linspace(0, 2 * np.pi, 100)
        data_y, data_x = bh.numpy.histogram(
            _data.phi.to_numpy(), bins=bins, range=(0, 2 * np.pi))
        x = (data_x[1:] + data_x[:-1]) / 2.0
        mc_rec_y, _ = bh.numpy.histogram(
            _mc_rec_data.phi.to_numpy(), bins=bins, range=(0, 2 * np.pi)
        )
        thrown_y, _ = bh.numpy.histogram(
            _thrown_data.phi.to_numpy(), bins=bins, range=(0, 2 * np.pi)
        )

        # Change 0's to 1 for division
        thrown_y = np.where(thrown_y == 0, 1, thrown_y)
        mc_rec_y = np.where(mc_rec_y == 0, 1, mc_rec_y)

        acceptance = np.nan_to_num(thrown_y / mc_rec_y)
        data_y = data_y / np.max(data_y)
        y = data_y * acceptance
        y = y / np.max(y)

        error_bar = np.ones_like(y)

        F = (mc_rec_y/thrown_y)
        error = np.sqrt(((thrown_y-mc_rec_y)*mc_rec_y) /
                        np.power(thrown_y, 3))/F
        error_bar = np.sqrt(
            np.power((y*error), 2) + np.power(stats.sem(y), 2))

        ax[a][b].set_ylim(bottom=0, top=np.max(y)*1.5)
        ax[a][b].errorbar(
            x,
            y,
            yerr=error_bar,
            marker=".",
            linestyle="",
            c="k",
            zorder=1,
        )

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

        _ax = ax[a][b].twinx()
        _ax.plot(phis, crossSections, c='r', linestyle='dotted')
        _ax.set_ylim(bottom=0, top=np.max(crossSections)*1.5)

        popt, pcov = curve_fit(func, x, y, maxfev=8000)

        ax[a][b].plot(xs, func(xs, *popt), c="#9467bd",
                      linewidth=2.0)
        ax[a][b].text(0.0, 1.2, cos_label)
    if not os.path.exists(f'{out_folder}/CosT'):
        os.makedirs(f'{out_folder}/CosT')

    plt.savefig(
        f"{out_folder}/CosT/W[{np.round(w.left,3)},{np.round(w.right,3)}]_Q2[{np.round(q2.left,3)},{np.round(q2.right,3)}]_{bins}_CosT.png", bbox_inches='tight'
    )


def draw_cos_plots(func, data, mc_rec_data, thrown_data, w, q2, cos_t_bins, out_folder):
    # with multiprocessing.Pool() as pool:
    #     inputs = []
    #     for bins in range(10, 11):
    #         inputs.append((func, data, mc_rec_data, thrown_data,
    #                        w, q2, cos_t_bins, out_folder, bins))
    #     pool.starmap(draw_cos_bin, inputs)
    draw_cos_bin(func, data, mc_rec_data, thrown_data,
                 w, q2, cos_t_bins, out_folder, 10)


def draw_xsec_plots(func, data, mc_rec_data, thrown_data, w, q2, cos_t, out_folder, bins):
    fig, ax = plt.subplots(2, 2, figsize=(12, 9))
    fig.suptitle(
        f"W={w},\t$Q^2$={q2},\tcos($\\theta$)={cos_t} \n{bins} bins in $\phi$"
    )
    xs = np.linspace(0, 2 * np.pi, 100)
    data_y, data_x = bh.numpy.histogram(
        data.phi, bins=bins, range=(0, 2 * np.pi))
    x = (data_x[1:] + data_x[:-1]) / 2.0
    mc_rec_y, _ = bh.numpy.histogram(
        mc_rec_data.phi, bins=bins, range=(0, 2 * np.pi)
    )
    thrown_y, _ = bh.numpy.histogram(
        thrown_data.phi, bins=bins, range=(0, 2 * np.pi)
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

    popt, pcov = curve_fit(func, x, y, maxfev=8000)
    # To compute one standard deviation errors on the parameters use
    # https://stackoverflow.com/questions/49130343/is-there-a-way-to-get-the-error-in-fitting-parameters-from-scipy-stats-norm-fit
    perr = np.sqrt(np.diag(pcov))
    ax[1][1].plot(xs, func(xs, *popt), c="#9467bd",
                  linewidth=2.0)

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


def draw_xsection(rec, mc_rec, thrown, func, out_folder, wbins, q2bins, thetabins):
    # executor = get_reusable_executor(max_workers=len(np.unique(rec.theta_bin)))
    total_num = (
        len(wbins)
        * len(q2bins)
    )

    pbar = tqdm(total=total_num)

    for w in wbins:
        for q2 in q2bins:
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
            draw_cos_plots(func, data, mc_rec_data,
                           thrown_data, w, q2, rec.theta_bin, out_folder)
            #################################################
            for cos_t in thetabins:
                rec_cut = (
                    (w == rec.w_bin) & (q2 == rec.q2_bin) & (
                        cos_t == rec.theta_bin)
                )
                mc_rec_cut = (
                    (w == mc_rec.w_bin)
                    & (q2 == mc_rec.q2_bin)
                    & (cos_t == mc_rec.theta_bin)
                )
                thrown_cut = (
                    (w == thrown.w_bin)
                    & (q2 == thrown.q2_bin)
                    & (cos_t == thrown.theta_bin)
                )

                data = rec[rec_cut]
                mc_rec_data = mc_rec[mc_rec_cut]
                thrown_data = thrown[thrown_cut]

                # results = executor.map(
                #     draw_plots,
                #     (func, data, mc_rec_data, thrown_data, w, q2, cos_t, out_folder),
                # )

                draw_plots(
                    func, data, mc_rec_data, thrown_data, w, q2, cos_t, out_folder
                )
            pbar.update(1)

    pbar.close()


def draw_kinematics(rec, w_bins, q2_bins, theta_bins, name="reconstructed"):
    rec = rec.dropna()

    fig, ax = plt.subplots(figsize=(12, 9))
    h = ax.hist2d(rec.w.to_numpy(), rec.q2.to_numpy(),
                  bins=200, range=[[np.min(w_bins), np.max(w_bins)], [np.min(q2_bins), np.max(q2_bins)]])

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
    ax2.set_xlabel(f"W [$\mathrm{{{{GeV}}}}$]", fontsize=30)
    ax2.tick_params(axis='both', which='major', labelsize=15)
    ax2.set_xlim(np.min(w_bins), np.max(w_bins))
    ax2.set_ylim(bottom=0)

    fig2.suptitle(f"W (n$\pi^+$)", fontsize=30)
    fig2.savefig(f"{out_folder}/kinematics/W_{name}.png", bbox_inches='tight')
    #########################################

    fig1, ax1 = plt.subplots(figsize=(12, 9))
    h = ax1.hist2d(rec.cos_theta.to_numpy(), rec.phi.to_numpy(), bins=100)
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
        "--e1f",
        action='store_true'
    )

    args = parser.parse_args()
    if args.e1f:
        ENERGY = 5.479

    pyarrow.set_cpu_count(8)

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
    mc_rec["cos_theta"] = np.cos(mc_rec.theta)

    mc_thrown = mc_thrown[(mc_thrown.w > 0)]
    mc_thrown["cos_theta"] = np.cos(mc_thrown.theta)

    start = time.time()
    rec = feather.read_feather(rec_data_file_path)
    stop = time.time()
    # print(f"\n\nread time rec: {stop - start}\n\n")
    # rec = rec[(rec.w > 0) & (rec.mm2 > 0.5) & (rec.mm2 < 1.5)]
    rec["cos_theta"] = np.cos(rec.theta).astype(np.float32)
    # print(f"===========================\nmc_rec:\n\n")
    # print(f"{mc_rec.info(verbose=True, memory_usage='deep')}")
    # print(f"\n\n===========================")
    # print(f"===========================\nmc_thrown:\n\n")
    # print(f"{mc_thrown.info(verbose=True, memory_usage='deep')}")
    # print(f"\n\n===========================")
    # print(f"===========================\nrec:\n\n")
    # print(f"{rec.info(verbose=True, memory_usage='deep')}")
    # print(f"\n\n===========================")

    sector_cuts = mm_cut(rec)

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
                     "phi", "helicty"]].copy(deep=True)
    mc_thrown = mc_thrown[["w", "q2", "mm2", "cos_theta", "phi", "helicty"]].copy(
        deep=True
    )
    rec = rec[["w", "q2", "mm2", "cos_theta",
               "phi", "helicty", "photon_flux"]].copy(deep=True)

    # Specifically put in bin edges
    # TODO ##################### BINS ######################
    w_bins = np.array([1.1, 1.125, 1.15, 1.175, 1.2, 1.225, 1.25, 1.275, 1.3,
                       1.325, 1.35, 1.375, 1.4, 1.425, 1.45, 1.475, 1.5, 1.525,
                       1.55, 1.575, 1.6, 1.625, 1.65, 1.675, 1.7, 1.725, 1.75,
                       1.775, 1.8])
    q2_bins = np.array([1.0, 1.4, 1.8, 2.6, 3.5])
    theta_bins = np.array([-1.0, -0.8, -0.6, -0.4, -0.2,
                           0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2])

    # TODO ##################### BINS ######################

    draw_kinematics(rec, w_bins, q2_bins, theta_bins)
    draw_kinematics(mc_rec, w_bins, q2_bins, theta_bins, "mc_rec")
    draw_kinematics(mc_thrown, w_bins, q2_bins, theta_bins, "thrown")

    mc_rec["w_bin"] = pd.cut(mc_rec["w"], bins=w_bins, include_lowest=True)
    mc_rec["q2_bin"] = pd.cut(mc_rec["q2"], bins=q2_bins, include_lowest=True)
    mc_rec["theta_bin"] = pd.cut(
        mc_rec["cos_theta"], bins=theta_bins, include_lowest=True
    )

    mc_thrown["w_bin"] = pd.cut(
        mc_thrown["w"], bins=w_bins, include_lowest=True)
    mc_thrown["q2_bin"] = pd.cut(
        mc_thrown["q2"], bins=q2_bins, include_lowest=True)
    mc_thrown["theta_bin"] = pd.cut(
        mc_thrown["cos_theta"], bins=theta_bins, include_lowest=True
    )

    rec["w_bin"] = pd.cut(rec["w"], bins=w_bins, include_lowest=True)
    rec["q2_bin"] = pd.cut(rec["q2"], bins=q2_bins, include_lowest=True)
    rec["theta_bin"] = pd.cut(
        rec["cos_theta"], bins=theta_bins, include_lowest=True)

    mc_rec.dropna(inplace=True)
    mc_thrown.dropna(inplace=True)
    rec.dropna(inplace=True)

    wbins = pd.unique(rec.w_bin)
    q2bins = pd.unique(rec.q2_bin)
    thetabins = pd.unique(rec.theta_bin)

    draw_xsection(rec, mc_rec, mc_thrown, model,
                  out_folder, wbins, q2bins, thetabins)

    stop = time.time()
    print(f"\n\nFull Running time: {stop - total_time}\n\n")
    print(
        f"\n\nFull Running time: {datetime.timedelta(seconds=(time.time() - total_time))}\n\n"
    )
