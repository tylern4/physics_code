#!/usr/bin/env python

import matplotlib  # noqa
matplotlib.use("agg")  # noqa
import warnings  # noqa
warnings.filterwarnings("ignore")  # noqa

from typing import Dict
import lmfit
# from loky import get_reusable_executor
import multiprocessing
import os
from maid_interface import maid_2007_Npi as maid
import datetime
import boost_histogram as bh
from pyarrow import csv, feather
import time
import argparse
from scipy.optimize import curve_fit
from scipy import stats
from scipy.special import erfc
from tqdm import tqdm
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from lmfit import Model, Parameters
from lmfit.models import *
from calc_xsections import *


ENERGY = 4.81726


def read_csv(file_name: str = "", data: bool = False):
    names = [
        "electron_sector",
        "w",
        "q2",
        "theta",
        "phi",
        "mm2",
        "cut_fid",
        "helicty",
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
        "cut_fid": "bool",
    }

    start = time.time()

    # Load file into pyTable before changing to pandas
    pyTable = csv.read_csv(
        file_name,
        read_options=csv.ReadOptions(
            use_threads=True, column_names=names),
        convert_options=csv.ConvertOptions(column_types=dtype),
    )
    df = pyTable.to_pandas(strings_to_categorical=True)

    if data:
        return df.dropna()

    mc_rec = df[df.type == "mc_rec"].dropna()
    thrown = df[df.type == "thrown"].dropna()
    del df

    stop = time.time()

    return (
        mc_rec,
        thrown,
    )


# def model(x, a, b, c):
#     """
#     a => sigma_l + sigma_t
#     b => epsilon*sigma_tt
#     c => Sqrt(2epsilon(1+epsilon))* sigma_lt
#     """
#     f = a + b * np.cos(2 * x) + c * np.cos(x)
#     return f


def A(M, B, C):
    if (C > 0 and np.abs(B) <= 4*C):
        return M**2 + B**2/(8*C) + C
    else:
        return M**2 + np.abs(B) - C


def model_new(x, M, b, c):
    """
    a => sigma_l + sigma_t
    b => epsilon*sigma_tt
    c => Sqrt(2epsilon(1+epsilon))* sigma_lt
    """
    f = A(M, b, c) + b * np.cos(2*x) + c * np.cos(x)
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


def virtual_photon(W: float, Q2: float, beam_energy: float) -> float:
    MASS_E = 0.000511
    target_mass = 0.93827203
    FS_ALPHA = 0.007297352570866302

    one = FS_ALPHA/(4 * np.pi)
    two = W/(beam_energy**2 * target_mass**2 * Q2)
    three = (W**2 - target_mass**2)

    beam_momentum = np.sqrt(beam_energy**2 - MASS_E**2)
    nu = (((W**2 + Q2) / target_mass) - target_mass) / 2  # Photon Energy
    scattered_energy = (beam_energy - nu)
    scattered_momentum = np.sqrt(scattered_energy**2 - MASS_E**2)
    theta = np.arccos((beam_energy * scattered_energy - Q2 / 2.0 - MASS_E**2) /
                      (beam_momentum * scattered_momentum))
    epsilon = 1/(1 + (2 * (1 + ((nu**2) / Q2)) * np.tan(theta / 2)**2))

    four = 1/(1 - epsilon)

    # This makes it look closer? Where am I off?
    flux = one * two * three * four * 10**3

    return flux


def hist_data(data, density=True, bins=10):
    data_y, data_x = bh.numpy.histogram(
        data.phi.to_numpy(), bins=bins, range=(0, 2 * np.pi), density=density, threads=4)
    x = (data_x[1:] + data_x[:-1]) / 2.0
    return data_y, x


def momentum_fn(energy, mass):
    return np.sqrt(energy**2 - mass**2)


def virtual_photon_energy_fn(target_mass, w, q2):
    return (((w**2 + q2) / target_mass) - target_mass) / 2


def costheta_e_fn(beam_energy, target_mass, w, q2):
    electron_mass = 5.109989433549345e-4
    nu = virtual_photon_energy_fn(target_mass, w, q2)
    scattered_energy = ((beam_energy) - (nu))
    beam_momentum = momentum_fn(beam_energy, electron_mass)
    scattered_momentum = momentum_fn(scattered_energy, electron_mass)
    return (beam_energy * scattered_energy - q2 / 2.0 - electron_mass**2) / (beam_momentum * scattered_momentum)


def virtual_photon_epsilon_fn(beam_energy, target_mass, w, q2):
    theta_e = np.arccos(costheta_e_fn(beam_energy, target_mass, w, q2))
    nu = virtual_photon_energy_fn(target_mass, w, q2)
    return np.power(1 + 2 * (1 + nu**2 / q2) * np.power(np.tan(theta_e / 2), 2), -1)


def virtual_photon_flux(w: float, q2: float, beam_energy: float = 4.81726, target_mass: float = 0.93827203) -> float:
    alpha = 0.007297352570866302
    epsilon = virtual_photon_epsilon_fn(beam_energy, target_mass, w, q2)
    return alpha / (4 * np.pi * q2) * w / (beam_energy**2 * target_mass**2) * (w**2 - target_mass**2) / (1 - epsilon)


def mm_cut(df: pd.DataFrame, sigma: int = 4, lmfit_fitter: bool = False) -> Dict:
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

        plt.errorbar(x, y, yerr=stats.sem(y), fmt=".", zorder=1)
        ax[a][b].errorbar(
            x, y, yerr=stats.sem(y), fmt=".", zorder=1)

        if lmfit_fitter:
            peak = PseudoVoigtModel(prefix="peak_")
            pars = peak.guess(y, x=x)
            background = GaussianModel(prefix="back_")
            pars.update(background.make_params())
            model = peak + background

            out = model.fit(y, pars, x=x)
            xs = np.linspace(0.3, 1.5, 1000)
            comps = out.eval_components(x=xs)
            out.params.pretty_print()
            ys = out.eval(params=out.params, x=xs)

            plt.plot(xs, comps['peak_'],
                     'r-', label='Peak Component')
            plt.plot(xs, comps['back_'],
                     'k--', label='Background Component')

            plt.plot(xs, ys, 'r-', linewidth=2.0, alpha=0.4,
                     label=f"Peak Center: {out.params['peak_center'].value:0.4f}")
            plt.axvline(out.params['peak_center']+sigma *
                        out.params['peak_fwhm'] / 2.355, c='r', alpha=0.4)
            plt.axvline(out.params['peak_center']-sigma *
                        out.params['peak_fwhm'] / 2.355, c='r', alpha=0.4)
            ax[a][b].plot(xs, ys, 'r-', linewidth=2.0,
                          alpha=0.6, label=f"Sector {sec}")
            ax[a][b].axvline(out.params['peak_center']+sigma *
                             out.params['peak_fwhm'] / 2.355, c='r', alpha=0.6)
            ax[a][b].axvline(out.params['peak_center']-sigma *
                             out.params['peak_fwhm'] / 2.355, c='r', alpha=0.6)
            data[sec] = (out.params['peak_center']-sigma*out.params['peak_fwhm'] / 2.355,
                         out.params['peak_center']+sigma*out.params['peak_fwhm'] / 2.355)
        else:
            popt_g, pcov_g = curve_fit(gauss, x, y, maxfev=8000)
            plt.plot(x, gauss(x, *popt_g), linewidth=2.0, alpha=0.4)

            p0 = [popt_g[0], popt_g[1], popt_g[2], 1.0, 1.0]
            popt, pcov = curve_fit(degauss, x, y, maxfev=8000)

            plt.plot(x, degauss(x, *popt), c="#9467bd",
                     linewidth=2.0, alpha=0.4)

            ax[a][b].plot(x, degauss(
                x, *popt), c="#9467bd", linewidth=3.0, label=f"Sector {sec}")

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
            data[sec] = (popt[1] - sigma * fwhm / 2.355,
                         popt[1] + sigma * fwhm / 2.355)

        ax[a][b].legend(loc='upper right')
        ax[a][b].set_xlabel(f"Mass [ $\mathrm{{{{GeV}}}}^2$]")

        if not os.path.exists(f'{out_folder}/cuts'):
            os.makedirs(f'{out_folder}/cuts')

        plt.xlabel(f"Mass [ $\mathrm{{{{GeV}}}}^2$]")
        plt.legend(loc='upper right')
        plt.title(
            f"Missing Mass Squared $e\left( p, \pi^{{{'+'}}} X \\right)$ in sector {sec}")

        plt.savefig(f"{out_folder}/cuts/MM2_cut_{sec}.png",
                    bbox_inches='tight')

    fig.suptitle(
        f"Missing Mass Squared $e\left( p, \pi^{{{'+'}}} X \\right)$", fontsize=20)
    fig.savefig(f"{out_folder}/cuts/MM2_cut_all.png",
                bbox_inches='tight')
    return data


def draw_cos_bin(data, mc_rec_data, thrown_data, w, q2, cos_t_bins, out_folder, bins, models_fits={"model": model_new}):
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
    # subplot_kw={'projection': 'polar'}

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

    phi_bins = np.linspace(0, 2 * np.pi, 200)
    thetabins = pd.unique(cos_t_bins)
    q2_left = 0.0
    for c, cos_t in enumerate(thetabins):
        a = which_plot[round(cos_t.left, 1)][0]
        b = which_plot[round(cos_t.left, 1)][1]
        cos_label = plot_label[round(cos_t.left, 1)]

        _w = (w.left + w.right) / 2.0
        _q2 = (q2.left + q2.right) / 2.0
        _cos_t = (cos_t.left + cos_t.right) / 2.0

        crossSections = []
        phis = []
        for phi in phi_bins:
            crossSections.append(
                maid(ENERGY, _w, _q2, _cos_t, np.degrees(phi)))
            phis.append(phi)

        crossSections = np.array(crossSections)
        phis = np.array(phis)

        _data = data[cos_t == data.theta_bin].copy()
        _mc_rec_data = mc_rec_data[cos_t == mc_rec_data.theta_bin].copy()
        _thrown_data = thrown_data[cos_t == thrown_data.theta_bin].copy()

        # flux = virtual_photon(w.left, q2.left, ENERGY)
        flux = virtual_photon_flux(w.left, q2.left)

        xs = np.linspace(0, 2 * np.pi, 100)
        # print(len(_data.cut_fid.to_numpy()), len(_data.phi.to_numpy()))
        cut_fids = {"With Fid cuts": True, "No fid cuts": False}
        for name, cuts in cut_fids.items():
            if cuts:
                # Histogram the data for plotting
                data_y, x = hist_data(
                    _data[_data.cut_fid], density=True, bins=bins)
                mc_rec_y, _ = hist_data(
                    _mc_rec_data[_mc_rec_data.cut_fid], density=True, bins=bins)
            else:
                # Histogram the data for plotting
                data_y, x = hist_data(data, density=True, bins=bins)
                mc_rec_y, _ = hist_data(_mc_rec_data, density=True, bins=bins)
                _ax = ax[a][b].twinx()
                _ax.plot(phis, crossSections, c='r', linestyle='dotted')
                _ax.set_ylim(bottom=0, top=np.max(crossSections)*1.5)
                _ax.text(0.0, 1.2*np.max(crossSections), cos_label)

            thrown_y, _ = hist_data(_thrown_data, density=True, bins=bins)

            # Change 0's to mean for division
            # thrown_y = thrown_y/np.max(thrown_y)
            # mc_rec_y = mc_rec_y/np.max(mc_rec_y)

            # Drop places with 0's
            cut = ~(mc_rec_y == 0) & ~(data_y == 0)
            x = x[cut]
            data_y = data_y[cut]
            thrown_y = thrown_y[cut]
            mc_rec_y = mc_rec_y[cut]

            acceptance = np.nan_to_num(thrown_y / mc_rec_y)

            try:
                data_y = data_y / np.max(data_y)
            except ValueError as e:
                print(e)
                continue

            y = (data_y * acceptance)  # * flux

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
                                  linewidth=2.0, c=ebar[0].get_color())
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
    for bins in range(10, 26, 2):
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

    sector_cuts = mm_cut(rec, sigma=10, lmfit_fitter=True)

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

    # mc_rec["w_bin"] = pd.cut(mc_rec["w"], bins=w_bins, include_lowest=False)
    # mc_rec["q2_bin"] = pd.cut(mc_rec["q2"], bins=q2_bins, include_lowest=False)
    # mc_rec["theta_bin"] = pd.cut(
    #     mc_rec["cos_theta"], bins=theta_bins, include_lowest=False
    # )

    # mc_thrown["w_bin"] = pd.cut(
    #     mc_thrown["w"], bins=w_bins, include_lowest=False)
    # mc_thrown["q2_bin"] = pd.cut(
    #     mc_thrown["q2"], bins=q2_bins, include_lowest=False)
    # mc_thrown["theta_bin"] = pd.cut(
    #     mc_thrown["cos_theta"], bins=theta_bins, include_lowest=False
    # )

    # rec["w_bin"] = pd.cut(rec["w"], bins=w_bins, include_lowest=False)
    # rec["q2_bin"] = pd.cut(rec["q2"], bins=q2_bins, include_lowest=False)
    # rec["theta_bin"] = pd.cut(
    #     rec["cos_theta"], bins=theta_bins, include_lowest=False)

    # mc_rec.dropna(inplace=True)
    # mc_thrown.dropna(inplace=True)
    # rec.dropna(inplace=True)

    print(f"===========================\nmc_rec:\n\n")
    print(f"{mc_rec.info(verbose=True, memory_usage='deep')}")
    print(f"\n\n===========================")
    print(f"===========================\nmc_thrown:\n\n")
    print(f"{mc_thrown.info(verbose=True, memory_usage='deep')}")
    print(f"\n\n===========================")
    print(f"===========================\nrec:\n\n")
    print(f"{rec.info(verbose=True, memory_usage='deep')}")
    print(f"\n\n===========================")

    # binning = dict()
    # binning["wbins"] = pd.Index.sort_values(pd.unique(rec.w_bin))
    # binning["q2bins"] = pd.Index.sort_values(pd.unique(rec.q2_bin))
    # binning["thetabins"] = pd.Index.sort_values(pd.unique(rec.theta_bin))

    # draw_xsection(rec, mc_rec, mc_thrown, model_new,
    #               out_folder, binning)

    stop = time.time()
    print(f"\n\nFull Running time: {stop - total_time}\n\n")
    print(
        f"\n\nFull Running time: {datetime.timedelta(seconds=(time.time() - total_time))}\n\n"
    )
