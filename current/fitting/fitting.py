#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import root_pandas

from iminuit import Minuit
from probfit import UnbinnedLH, BinnedChi2, BinnedLH, Extended, AddPdf, AddPdfNorm, SimultaneousFit
from probfit.pdf import crystalball, rtv_breitwigner, linear, poly2, poly3, gaussian
from probfit.plotting import draw_pdf, draw_normed_pdf, draw_compare_hist

def __fit1(binned_func, df, bound, bins, sector):
    print(sector)
    plt.figure(figsize=(16, 9))
    temp = df[df.sector == sector]
    fit_data = np.array(temp.MM, dtype=np.float64)
    model = AddPdf(gaussian, rtv_breitwigner, poly3)
    binned_function = binned_func(model, fit_data, bins=bins, bound=bound)
    minuit = Minuit(binned_function, mean=1.0, sigma=0.5, error_mean=0.1, error_sigma=0.1,
                                        m=0.94, gamma=0.03, error_m=0.1, error_gamma=0.001,
                                        a=0.5, b=0.2, c=0.1, d=0.1,
                                        error_a=0.01, error_b=0.01, error_c=0.01, error_d=0.01)
    minuit.migrad()
    binned_function.draw(minuit)
    plt.savefig(f"mm_fit_{binned_func.__name__}_{sector}.png")
    print(sector)
    return True

def __fit(binned_func, df, bound, bins, sector):
    print(sector)
    plt.figure(figsize=(16, 9))
    temp = df[df.sector == sector]
    fit_data = np.array(temp.MM, dtype=np.float64)
    model = AddPdfNorm(rtv_breitwigner)
    binned_function = binned_func(model, fit_data, bins=bins, bound=bound)
    minuit = Minuit(binned_function, m=0.94, gamma=0.03, error_m=0.1, error_gamma=0.001)
    minuit.migrad()
    try:
        minuit.minos()
    except Exception as e:
        print(e)
    try:
        minuit.hesse()
    except Exception as e:
        print(e)
    minuit.migrad()
    binned_function.draw(minuit)
    plt.savefig(f"mm_fit_{binned_func.__name__}_{sector}.png")
    print(sector)
    return True

def fit_mm_secBysec(binned_func, df, bound=(0.8, 1.1), bins=1500):
    #binned_func -> BinnedChi2, BinnedLH
    args = [__fit(binned_func, df, bound, bins,i) for i in range(1, 7)]
    return True


if __name__ == '__main__':
    bound=(0.8, 1.1)
    df = root_pandas.read_root("~/Data/e1d/ntuple.root", "ntuple")
    fit_mm_secBysec(BinnedLH, df, bound=bound, bins=200)
