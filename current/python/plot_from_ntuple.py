#!/usr/bin/env python

import pandas as pd
import root_pandas as root_pd
import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy.optimize import curve_fit
from scipy.stats import norm


def gaus(x, a, mu, sigma):
    return a * norm.pdf(x, mu, sigma)


def plot_elastic(data):
    data = data[data.type == 33]
    fig, axs = plt.subplots(2, 3, figsize=(16, 9))
    for i in range(0, 6):
        tmp = data[data.sector == i+1]
        axs[i % 2, i % 3].set_title(f"sector {i+1}")
        ydata, xdata, _ = axs[i % 2, i % 3].hist(
            tmp.W, bins=250, range=(0.5, 2.0), density=True)

        popt, pcov = curve_fit(gaus, xdata[1:], ydata, p0=(0.1, 0.94, 0.001))
        axs[i % 2, i % 3].plot(xdata, gaus(xdata, *popt), 'r--',
                               label=f'fit: {popt[1]:0.3f}, {popt[2]:0.3f}')

        axs[i % 2, i % 3].legend()

    plt.show()


def plot_mm(data):
    data = data[data.type == 0]
    fig, axs = plt.subplots(2, 3, figsize=(16, 9))
    for i in range(0, 6):
        tmp = data[data.sector == i+1]
        axs[i % 2, i % 3].set_title(f"sector {i+1}")
        ydata, xdata, _ = axs[i % 2, i % 3].hist(
            tmp.MM, bins=250, range=(0.5, 2.0), density=True)

        popt, pcov = curve_fit(gaus, xdata[1:], ydata, p0=(0.1, 0.94, 0.001))
        axs[i % 2, i % 3].plot(xdata, gaus(xdata, *popt), 'r--',
                               label=f'fit: {popt[1]:0.3f}, {popt[2]:0.3f}')

        axs[i % 2, i % 3].legend()

    plt.show()


if __name__ == "__main__":
    if len(sys.argv) > 1:
        data = root_pd.read_root(sys.argv[1])
    else:
        data = root_pd.read_root("/Users/tylern/Data/ntuple/ntuple_1.root")

    # plot_elastic(data)
    plot_mm(data)
