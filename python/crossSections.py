#!/usr/bin/env python

import matplotlib  # noqa
matplotlib.use("agg")  # noqa
import warnings  # noqa
warnings.filterwarnings("ignore")  # noqa

import matplotlib.pyplot as plt
import inspect
from lmfit import Model
from lmfit.models import *
from maid_interface import maid_2007_Npi as maid
import argparse
import numpy as np
import pandas as pd
import pyarrow as pa
from pyarrow import csv
import time
from tqdm import tqdm
from scipy.special import erfc
from scipy import stats
from cycler import cycler
import boost_histogram as bh
import os


def degauss(x, A, mu, sigma, lambda1, lambda2):
    mu1 = sigma * sigma * lambda1 + x - mu
    mu2 = -sigma * sigma * lambda2 + x - mu
    ret = A * 0.5 / (1.0 / lambda1 + 1.0 / lambda2) * \
        (np.exp(0.5 * np.power(sigma * lambda1, 2) + lambda1 * (x - mu)) * erfc(mu1 / (sigma * np.sqrt(2.0)))
         + np.exp(0.5 * np.power(sigma * lambda2, 2) - lambda2 * (x - mu)) * erfc(-mu2 / (sigma * np.sqrt(2.0))))

    return ret


def gauss(x, A, mu, sig):
    ret = np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))
    return A*ret


def peak(x, c):
    return np.exp(-np.power(x - c, 2) / 16.0)


def lin_interp(x, y, i, half):
    return x[i] + (x[i+1] - x[i]) * ((half - y[i]) / (y[i+1] - y[i]))


def half_max_x(x, y):
    half = np.max(y)/2.0
    signs = np.sign(np.add(y, -half))
    zero_crossings = (signs[0:-2] != signs[1:-1])
    zero_crossings_i = np.where(zero_crossings)[0]
    return [lin_interp(x, y, zero_crossings_i[0], half),
            lin_interp(x, y, zero_crossings_i[1], half)]


def mm_cut(df: pd.DataFrame, sigma: int = 4):
    data = {}

    for sec in range(1, 7):
        y, x = bh.numpy.histogram(
            df[df.electron_sector == sec].mm2, bins=500, density=True
        )
        x = (x[1:] + x[:-1]) / 2

        peak = PseudoVoigtModel(prefix="peak_")
        pars = peak.guess(y, x=x)
        background = GaussianModel(prefix="back_")
        pars.update(background.make_params())
        model = peak * background
        out = model.fit(y, pars, x=x)
        xs = np.linspace(0.3, 1.5, 1000)
        data[sec] = (out.params['peak_center']-sigma*out.params['peak_fwhm'] / 2.355,
                     out.params['peak_center']+sigma*out.params['peak_fwhm'] / 2.355)

    return data


def cut_for_MM(rec, mc_rec):
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

    return rec, mc_rec


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
        "cut_fid": "bool"
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


def hist_data(data, density=True):
    data_y, data_x = bh.numpy.histogram(
        data.phi.to_numpy(), bins=10, range=(0, 2 * np.pi), density=True, threads=4)
    x = (data_x[1:] + data_x[:-1]) / 2.0
    return data_y, x


def prep_for_ana(dataframe, w_bins, q2_bins, theta_bins):
    dataframe.dropna(inplace=True)
    dataframe["cos_theta"] = np.cos(dataframe.theta).astype(np.float32)
    dataframe["w_bin"] = pd.cut(
        dataframe["w"], bins=w_bins, include_lowest=False)
    dataframe["q2_bin"] = pd.cut(
        dataframe["q2"], bins=q2_bins, include_lowest=False)
    dataframe["theta_bin"] = pd.cut(
        dataframe["cos_theta"], bins=theta_bins, include_lowest=False)
    dataframe.dropna(inplace=True)

    return dataframe


def make_cuts(dataframe, w, q2, theta):
    return (dataframe.w_bin == w) & (dataframe.q2_bin == q2) & (dataframe.theta_bin == theta)


def get_maid_values(xs, w, q2, theta):
    ENERGY = 4.81726
    crossSections = []
    for phi in xs:
        crossSections.append(maid(ENERGY, w, q2, theta, np.degrees(phi)))

    return np.array(crossSections)


def get_error_bars(y, mc_rec_y, thrown_y):
    F = (mc_rec_y/thrown_y)
    error = np.sqrt(((thrown_y-mc_rec_y)*mc_rec_y) / np.power(thrown_y, 3))/F
    error_bar = np.sqrt(np.power((y*error), 2) + np.power(stats.sem(y), 2))

    return error_bar


def plot_maid_model(ax, w, q2, theta, xs):
    # Get bin centers
    _w = (w.left + w.right) / 2.0
    _q2 = (q2.left + q2.right) / 2.0
    _theta = (theta.left + theta.right) / 2.0
    # Get the cross section values from maid
    crossSections = get_maid_values(xs, _w, _q2, _theta)
    _ax = ax.twinx()
    _ax.plot(xs, crossSections, c='r', linestyle='dotted')
    _ax.set_ylim(bottom=0, top=np.max(crossSections)*1.5)

    return _ax


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


def fit_model(ax, func, x, y, xs, color, name):
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
    except ValueError as e:
        print(e)
        return None
    except TypeError as e:
        print(e)
        return None

    # Plot the fitted model with output parameters and same x's as model
    ax.plot(xs, out.eval(params=out.params, x=xs),
            linewidth=2.0, c=color, label=f'{name}')

    # Get uncertinty and plot between 2 sigmas
    dely = out.eval_uncertainty(sigma=2, x=xs)
    ax.fill_between(xs, out.eval(params=out.params, x=xs)-dely,
                    out.eval(
        params=out.params, x=xs)+dely,
        color=color, alpha=0.1)

    return out


def main(rec, mc_rec, mc_thrown, binning, out_folder="plots"):
    if not os.path.exists(f'{out_folder}/crossSections'):
        os.makedirs(f'{out_folder}/crossSections')
    # Make a set of values from 0 to 2Pi for plotting
    xs = np.linspace(0, 2 * np.pi, 250)
    for w in tqdm(binning["wbins"]):
        for q2 in binning["q2bins"]:
            for theta in binning["thetabins"]:
                fig = plt.figure(figsize=(12, 9), constrained_layout=True)
                gs = fig.add_gridspec(2, 1, height_ratios=[3, 1])
                ax1 = fig.add_subplot(gs[0])
                ax2 = fig.add_subplot(gs[1])

                plot_maid_model(ax1, w, q2, theta, xs)

                # Cut data/mc for the w/q2/theta bin we're in
                data = rec[make_cuts(rec, w, q2, theta)].copy()
                data_mc = mc_rec[make_cuts(mc_rec, w, q2, theta)].copy()
                thrown = mc_thrown[make_cuts(mc_thrown, w, q2, theta)].copy()

                cut_fids = {"fid cuts": True, "no fid cuts": False}
                for name, cut in cut_fids.items():
                    if cut:
                        # Histogram the data for plotting
                        data_y, x = hist_data(data[data.cut_fid], density=True)
                        mc_rec_y, _ = hist_data(
                            data_mc[data_mc.cut_fid], density=True)
                        thrown_y, _ = hist_data(
                            thrown[thrown.cut_fid], density=True)
                    else:
                        # Histogram the data for plotting
                        data_y, x = hist_data(data, density=True)
                        mc_rec_y, _ = hist_data(data_mc, density=True)
                        thrown_y, _ = hist_data(thrown, density=True)

                    cut = ~(mc_rec_y == 0)
                    x = x[cut]
                    data_y = data_y[cut]
                    thrown_y = thrown_y[cut]
                    mc_rec_y = mc_rec_y[cut]

                    # Calculate acceptance and correct data
                    acceptance = np.nan_to_num(thrown_y / mc_rec_y)
                    ax2.errorbar(x, 1/acceptance, marker=".",
                                 linestyle="", zorder=1, label=f"{name}")

                    y = (data_y * acceptance)
                    # Calc errorbars
                    error_bar = get_error_bars(y, mc_rec_y, thrown_y)
                    ebar = ax1.errorbar(
                        x, y, yerr=error_bar, marker=".", linestyle="", zorder=1, label=f"{name}")
                    out = fit_model(ax1, model_new, x, y, xs,
                                    ebar[0].get_color(), name)
                    ax1.legend()

                    top = np.max(y)*1.5
                    if np.isnan(top):
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
    args = parser.parse_args()

    # Start to main

    # Specifically put in bin edges
    # TODO ##################### BINS ######################
    w_bins = np.array([1.1, 1.125, 1.15, 1.175, 1.2, 1.225, 1.25, 1.275, 1.3,
                       1.325, 1.35, 1.375, 1.4, 1.425, 1.45, 1.475, 1.5, 1.525,
                       1.55, 1.575, 1.6, 1.625, 1.65, 1.675, 1.7, 1.725, 1.75,
                       1.775, 1.8])
    # w_bins = np.array([1.1, 1.2, 1.3, 1.4, 1.5, 1.6,  1.7,  1.8])
    # q2_bins = np.array([1.0, 1.4, 1.8, 2.6, 3.5])
    q2_bins = np.array([1.2, 1.6, 1.8, 2.2, 2.6, 3.5])
    theta_bins = np.array([-1.0, -0.8, -0.6, -0.4, -0.2,
                           0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2])
    # TODO ##################### BINS ######################

    print("Start setup")
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
    print("Done setup")

    main(rec, mc_rec, mc_thrown, binning, out_folder=args.out_folder)
