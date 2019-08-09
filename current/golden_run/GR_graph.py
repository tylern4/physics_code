#!/usr/bin/env python
from __future__ import print_function
import matplotlib.pyplot as plt
from scipy.stats import norm
import matplotlib.mlab as mlab
import argparse
import sys
from scipy import optimize
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('agg')

SIGMA = 12


def gaussian(x, amplitude, mean, stddev):
    return amplitude * np.exp(-((x - mean) / 4 / stddev)**2)


try:
    from termcolor import cprint
    def print_green(x): return cprint(x, 'green', attrs=['bold'])
    def print_red(x): return cprint(x, 'red', attrs=['bold'])
    def print_blue(x): return cprint(x, 'blue', attrs=['bold'])
    def print_white(x): return cprint(x, 'white', attrs=['bold'])
    def print_yellow(x): return cprint(x, 'yellow', attrs=['bold'])
except ImportError:
    def print_green(x): return print(x)
    def print_red(x): return print(x)
    def print_blue(x): return print(x)
    def print_white(x): return print(x)
    def print_yellow(x): return print(x)


def main():
    parser = argparse.ArgumentParser(
        description="Plot data from golden run output csv:")
    parser.add_argument('input', type=str, help="Input CSV")
    parser.add_argument(
        'output',
        type=str,
        nargs='?',
        help="Output folder for final root files",
        default='.')

    args = parser.parse_args()
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()
    if args.input[-3:] != 'csv':
        print_white(args.input[-3])
        print_red("Error: Input must be CSV file")
        parser.print_help()
        parser.exit()
    if args.output[-1] != '/':
        args.output = args.output + '/'

    golden_df = pd.read_csv(args.input)
    golden_df['ratio'] = golden_df['num_of_events'] / golden_df['total_q']

    fig = plt.figure(
        num=None, figsize=(16, 9), dpi=200, facecolor='w', edgecolor='k')
    bin_heights, bin_borders, _ = plt.hist(
        golden_df['ratio'], bins=100, histtype='stepfilled', alpha=0.9, range=(45000, 75000))
    popt, _ = optimize.curve_fit(
        gaussian, bin_borders[:-1], bin_heights, p0=[10.0, 60000.0, 200.0])
    plt.plot(bin_borders[:-1],
             gaussian(bin_borders[:-1], *popt), label=f"Mu={popt[1]:0.2f} Sigma={popt[2]:0.2f}")
    plt.axvline(x=popt[1])
    plt.axvline(x=popt[1]+SIGMA*popt[2])
    plt.axvline(x=popt[1]-SIGMA*popt[2])
    plt.legend()

    plt.xlim([40000, 80000])
    plt.xlabel('Ratio')
    plt.savefig('ratio_hist.pdf')

    fig = plt.figure(
        num=None, figsize=(16, 9), dpi=200, facecolor='w', edgecolor='k')
    plt.scatter(golden_df['run_num'], golden_df['ratio'], alpha=0.5)
    plt.axhline(y=popt[1])
    plt.axhline(y=popt[1]+SIGMA*popt[2])
    plt.axhline(y=popt[1]-SIGMA*popt[2])
    plt.xlabel('Run Number')
    plt.ylabel('Ratio (events/total Q)')
    plt.ylim([40000, 80000])
    plt.savefig('golden_run.pdf')

    keep = golden_df[(golden_df['ratio'] > popt[1] - SIGMA *
                      popt[2]) & (golden_df['ratio'] < popt[1]+SIGMA*popt[2])]
    keep = keep[['run_num', 'file_num']]

    # print(len(keep))
    # for index, row in keep.iterrows():
    #    print(
    #        "rsync -rav --partial --progress ../h10_r" + str(row[0]) + "_" +
    #       str(row[1]).zfill(2) + ".root",
    #        end=" .; \n")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print_red("\n\nExiting")
        sys.exit()
