#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import root_pandas as rp
import pandas as pd
from tqdm import tqdm
from scipy.special import erfc
import argparse
from scipy import stats
from scipy.optimize import curve_fit
import argparse

import warnings

warnings.filterwarnings("ignore")


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


def mm_cut(df):
    NSIGMA = 3
    data = {}
    for sec in range(1, 7):
        plt.figure(figsize=(12, 9))
        y, x = np.histogram(df[df.electron_sector == sec].mm2, bins=500, density=True)
        x = (x[1:] + x[:-1]) / 2
        popt_g, pcov_g = curve_fit(gauss, x, y, maxfev=8000)
        plt.plot(x, gauss(x, *popt_g), linewidth=2.0)
        plt.errorbar(x, y, yerr=np.sqrt(y.shape[0]) / y.shape[0], fmt=".", zorder=1)

        plt.axvline(popt_g[1] + NSIGMA * popt_g[2])
        plt.axvline(popt_g[1] - NSIGMA * popt_g[2])

        p0 = [popt_g[0], popt_g[1], popt_g[2], 1.0, 1.0]
        popt, pcov = curve_fit(degauss, x, y, maxfev=8000)

        plt.plot(x, degauss(x, *popt), c="#9467bd", linewidth=2.0)

        # find the FWHM
        xs = np.linspace(0.7, 1.5, 100000)
        hmx = half_max_x(xs, degauss(xs, *popt))
        fwhm = hmx[1] - hmx[0]
        plt.axvline(popt[1] + NSIGMA * fwhm / 2.355, c="#9467bd")
        plt.axvline(popt[1] - NSIGMA * fwhm / 2.355, c="#9467bd")

        data[sec] = (popt_g[1] + NSIGMA * popt_g[2], popt_g[1] - NSIGMA * popt_g[2])

        # print("{", end="")
        # for x in popt_g:
        #     print(f" {x:.20f},", end="")
        # print("}")

    return data


def draw_xsection(rec, mc_rec, thrown, func):
    xs = np.linspace(0, 2 * np.pi, 100)
    total_num = (
        len(np.unique(rec.w_bin))
        * len(np.unique(rec.q2_bin))
        * len(np.unique(rec.theta_bin))
    )
    pbar = tqdm(total=total_num)
    for w in np.unique(rec.w_bin):
        for q2 in np.unique(rec.q2_bin):
            for cos_t in np.unique(rec.theta_bin):
                pbar.update(1)
                rec_cut = (
                    (w == rec.w_bin) & (q2 == rec.q2_bin) & (cos_t == rec.theta_bin)
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
                fig, ax = plt.subplots(2, 2, figsize=(12, 9))
                fig.suptitle(f"W={w},\t$Q^2$={q2},\tcos($\Theta$)={cos_t}")
                for bins in range(10, 11):
                    data_y, data_x = np.histogram(
                        data.phi, bins=bins, range=(0, 2 * np.pi)
                    )
                    x = (data_x[1:] + data_x[:-1]) / 2.0
                    mc_rec_y, _ = np.histogram(
                        mc_rec_data.phi, bins=bins, range=(0, 2 * np.pi)
                    )
                    thrown_y, _ = np.histogram(
                        thrown_data.phi, bins=bins, range=(0, 2 * np.pi)
                    )

                    ax[0][0].errorbar(
                        x,
                        thrown_y,
                        marker=".",
                        yerr=stats.sem(thrown_y),
                        c="r",
                        linestyle="",
                        label="thrown",
                    )
                    ax[0][0].errorbar(
                        x,
                        mc_rec_y,
                        marker=".",
                        yerr=stats.sem(mc_rec_y),
                        c="orange",
                        linestyle="",
                        label="mc_rec",
                    )
                    ax[0][1].errorbar(
                        x,
                        data_y,
                        yerr=stats.sem(data_y),
                        marker=".",
                        linestyle="",
                        label="data",
                    )

                    acceptance = thrown_y / mc_rec_y
                    ax[1][0].errorbar(
                        x,
                        acceptance,
                        yerr=stats.sem(acceptance),
                        marker=".",
                        c="g",
                        linestyle="",
                        label="acceptance",
                    )

                    y = np.nan_to_num(data_y * (acceptance))

                    popt, pcov = curve_fit(func, x, y, maxfev=8000)
                    ax[1][1].errorbar(
                        x,
                        y,
                        yerr=stats.sem(y),
                        marker=".",
                        linestyle="",
                        c="k",
                        zorder=1,
                        label="corrected",
                    )

                    plt.plot(xs, func(xs, *popt), c="#9467bd", linewidth=2.0)
                fig.legend()
                plt.savefig(
                    f"plots/W[{w.left},{w.right}]_Q2[{q2.left},{q2.right}]_cos(theta)[{cos_t.left},{cos_t.right}].png"
                )
    pbar.close()


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

    args = parser.parse_args()

    mc_data_file_path = args.mc_data_file_path
    rec_data_file_path = args.rec_data_file_path

    mc_df = pd.read_csv(mc_data_file_path, index_col=False, memory_map=True)
    mc_df = mc_df[(mc_df.w > 0) & (mc_df.mm2 > 0.5) & (mc_df.mm2 < 1.5)]
    mc_df["cos_theta"] = np.cos(mc_df.theta)

    mc_rec = mc_df[mc_df.type == "mc_rec"]
    mc_thrown = mc_df[mc_df.type == "thrown"]

    mc_rec = mc_rec.merge(
        mc_thrown, left_on="hash", right_on="hash", suffixes=("", "_thrown")
    )
    mc_rec.drop(
        ["type", "hash", "type_thrown", "electron_sector_thrown", "helicty_thrown"],
        axis=1,
        inplace=True,
    )
    mc_thrown.drop(["type", "hash"], axis=1, inplace=True)

    rec = pd.read_csv(rec_data_file_path, index_col=False, memory_map=True)
    rec = rec[(rec.w > 0) & (rec.mm2 > 0.5) & (rec.mm2 < 1.5)]
    rec["cos_theta"] = np.cos(rec.theta)

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

    mc_rec = mc_rec[["w", "q2", "mm2", "cos_theta", "phi", "helicty"]].copy(deep=True)
    mc_thrown = mc_thrown[["w", "q2", "mm2", "cos_theta", "phi", "helicty"]].copy(
        deep=True
    )
    rec = rec[["w", "q2", "mm2", "cos_theta", "phi", "helicty"]].copy(deep=True)

    rec.head()

    w_bins = np.arange(1.2, 1.8, 0.025)
    q2_bins = np.arange(1.0, 3.5, 0.5)
    theta_bins = np.arange(-1.0, 1.0, 0.25)

    mc_rec["w_bin"] = pd.cut(mc_rec["w"], bins=w_bins, include_lowest=True)
    mc_rec["q2_bin"] = pd.cut(mc_rec["q2"], bins=q2_bins, include_lowest=True)
    mc_rec["theta_bin"] = pd.cut(
        mc_rec["cos_theta"], bins=theta_bins, include_lowest=True
    )

    mc_thrown["w_bin"] = pd.cut(mc_thrown["w"], bins=w_bins, include_lowest=True)
    mc_thrown["q2_bin"] = pd.cut(mc_thrown["q2"], bins=q2_bins, include_lowest=True)
    mc_thrown["theta_bin"] = pd.cut(
        mc_thrown["cos_theta"], bins=theta_bins, include_lowest=True
    )

    rec["w_bin"] = pd.cut(rec["w"], bins=w_bins, include_lowest=True)
    rec["q2_bin"] = pd.cut(rec["q2"], bins=q2_bins, include_lowest=True)
    rec["theta_bin"] = pd.cut(rec["cos_theta"], bins=theta_bins, include_lowest=True)

    mc_rec.dropna(inplace=True)
    mc_thrown.dropna(inplace=True)
    rec.dropna(inplace=True)

    draw_xsection(rec, mc_rec, mc_thrown, model)

