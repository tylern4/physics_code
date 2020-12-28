#!/usr/bin/env python

import warnings  # noqa
warnings.filterwarnings("ignore")  # noqa

from typing import Dict
from multiprocessing import Pool
import os
import datetime
import boost_histogram as bh
from pyarrow import csv, feather
import time
import argparse
from tqdm import tqdm
import pandas as pd
import numpy as np
from lmfit.models import *


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

    mc_rec = df[df.type == "mc_rec"]
    thrown = df[df.type == "thrown"]
    del df

    stop = time.time()

    return (
        mc_rec,
        thrown,
    )


def mm_cut(df: pd.DataFrame, sigma: int = 4) -> Dict:
    data = {}

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

        data[sec] = (out.params['peak_center']-sigma*out.params['peak_fwhm'] / 2.355,
                     out.params['peak_center']+sigma*out.params['peak_fwhm'] / 2.355)

    return data


def get_yeilds(data, mc_rec_data, thrown_data, w, q2, cos_t, bins=10):
    data_y, data_x = bh.numpy.histogram(
        data.phi, bins=bins, range=(0, 2 * np.pi), threads=4)
    x = (data_x[1:] + data_x[:-1]) / 2.0
    mc_rec_y, _ = bh.numpy.histogram(
        mc_rec_data.phi, bins=bins, range=(0, 2 * np.pi), threads=4
    )
    thrown_y, _ = bh.numpy.histogram(
        thrown_data.phi, bins=bins, range=(0, 2 * np.pi), threads=4
    )

    kin_bin = {"w": w, "q2": q2, "cos_t": cos_t}

    rec_keys = [f"data_{_x}" for _x in range(bins)]
    rec_data = dict(zip(rec_keys, data_y))

    mc_rec_keys = [f"mc_rec_{_x}" for _x in range(bins)]
    mc_rec_data = dict(zip(mc_rec_keys, mc_rec_y))

    thrown_keys = [f"thrown_{_x}" for _x in range(bins)]
    thrown_data = dict(zip(thrown_keys, thrown_y))

    return {**kin_bin, **rec_data, **mc_rec_data, **thrown_data}


def yeilds(rec, mc_rec, thrown, binning):
    total_num = (
        len(binning["wbins"])
        * len(binning["q2bins"])
        * len(binning["thetabins"])
    )

    args = []
    total_data = []
    pbar = tqdm(total=total_num)
    for w in binning["wbins"]:
        for q2 in binning["q2bins"]:
            for cos_t in binning["thetabins"]:
                pbar.update(1)
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

                data = rec[rec_cut].copy()
                mc_rec_data = mc_rec[mc_rec_cut].copy()
                thrown_data = thrown[thrown_cut].copy()

                total_data.append(get_yeilds(
                    data, mc_rec_data, thrown_data, w, q2, cos_t))

    pbar.close()

    # with Pool(8) as p:
    #     total_data = p.starmap(get_yeilds, args)

    return pd.DataFrame(total_data)


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
        "--e1f",
        action='store_true'
    )
    parser.add_argument(
        "--folder",
        dest="out_folder",
        type=str,
        help="Folder for yeilds",
        required=False,
        default="yeilds",
    )

    args = parser.parse_args()
    if args.e1f:
        ENERGY = 5.479

    total_time = time.time()
    mc_data_file_path = args.mc_data_file_path
    rec_data_file_path = args.rec_data_file_path
    out_folder = args.out_folder

    mc_rec, mc_thrown = read_csv(mc_data_file_path)
    mc_rec["cos_theta"] = np.cos(mc_rec.theta)

    mc_thrown = mc_thrown[(mc_thrown.w > 0)]
    mc_thrown["cos_theta"] = np.cos(mc_thrown.theta)

    rec = feather.read_feather(rec_data_file_path)

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
                     "phi", "helicty", "electron_sector"]].copy(deep=True)
    mc_thrown = mc_thrown[["w", "q2", "mm2", "cos_theta", "phi", "helicty", "electron_sector"]].copy(
        deep=True
    )
    rec = rec[["w", "q2", "mm2", "cos_theta",
               "phi", "helicty",  "electron_sector"]].copy(deep=True)

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

    binning = dict()
    binning["wbins"] = pd.unique(rec.w_bin)
    binning["q2bins"] = pd.unique(rec.q2_bin)
    binning["thetabins"] = pd.unique(rec.theta_bin)

    out = yeilds(rec, mc_rec, mc_thrown, binning)
    if not os.path.exists(f'{out_folder}'):
        os.makedirs(f'{out_folder}')

    out.to_csv(f"{out_folder}/yeilds.csv")
    out.to_pickle(f"{out_folder}/yeilds.pkl")

    stop = time.time()
    print(f"\n\nFull Running time: {stop - total_time}\n\n")
    print(
        f"\n\nFull Running time: {datetime.timedelta(seconds=(time.time() - total_time))}\n\n"
    )
