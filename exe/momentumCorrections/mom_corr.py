#!/usr/bin/env python

import pyarrow  # noqa
pyarrow.set_cpu_count(16)  # noqa
# import matplotlib  # noqa
# matplotlib.use("agg")  # noqa
import warnings  # noqa
warnings.filterwarnings('ignore')  # noqa

from lmfit.models import *
from scipy.optimize import curve_fit
from scipy import stats
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyarrow import csv
from lmfit.models import GaussianModel, PolynomialModel, BreitWignerModel, LognormalModel, Model
from multiprocessing import Pool
import sys
from typing import Dict
from matplotlib.pyplot import legend


MASS_P = 0.93827208816
E0 = 4.81726
MASS_E = 0.00051099895


def clock_to_time(clock: float) -> str:
    millis = clock * 1000
    hours = int((millis/(1000*60*60)) % 24)
    seconds = int((millis/1000) % 60)
    minutes = int((millis/(1000*60)) % 60)
    hours = int((millis/(1000*60*60)) % 24)
    if hours > 1:
        return f"{hours}:{minutes}:{seconds}"
    elif minutes > 1:
        return f"{minutes} min {seconds} sec"
    elif seconds > 1:
        return f"{seconds} sec"
    else:
        return f"{millis} ms"


def timeit(method):
    def timed(*args, **kw):
        ts = time.time()
        result = method(*args, **kw)
        te = time.time()

        print(f'{method.__name__ }  {clock_to_time(te - ts)}')
        return result
    return timed


def vec(method):
    return np.vectorize(method, cache=True, otypes=["float16"])


@vec
def calc_W(e_p, e_theta, e_phi):
    px = e_p*np.sin(e_theta)*np.cos(e_phi)
    py = e_p*np.sin(e_theta)*np.sin(e_phi)
    pz = e_p*np.cos(e_theta)

    e_beam = np.array([0, 0, np.sqrt(E0**2-MASS_E**2), E0])
    e_prime = np.array([px, py, pz, np.linalg.norm([px, py, pz, MASS_E])])
    p_targ = np.array([0, 0, 0.0, MASS_P])
    temp = e_beam - e_prime + p_targ
    temp2 = temp[0]*temp[0]+temp[1]*temp[1]+temp[2]*temp[2]-temp[3]*temp[3]
    return np.sqrt(-temp2)


@vec
def q2_calc(e_p, e_theta, e_phi):
    px = e_p*np.sin(e_theta)*np.cos(e_phi)
    py = e_p*np.sin(e_theta)*np.sin(e_phi)
    pz = e_p*np.cos(e_theta)

    e_beam = np.array([0, 0, np.sqrt(E0**2-MASS_E**2), E0])
    e_prime = np.array([px, py, pz, np.linalg.norm([px, py, pz, MASS_E])])

    temp = e_beam - e_prime
    temp2 = temp[0]*temp[0]+temp[1]*temp[1]+temp[2]*temp[2]-temp[3]*temp[3]

    return temp2


@vec
def center_phi(phi, sec):
    sector = {
        1: 90,
        2: 30,
        3: -30,
        4: -90,
        5: -150,
        6: 150
    }

    return phi-sector[sec]


@vec
def Theta_e_calc(theta_p):
    return 2 * np.arctan(MASS_P/((E0+MASS_P)*np.tan(theta_p)))


def FitFunc(e_phi, theta_e,
            alpha_A, beta_A, gamma_A,
            alpha_B, beta_B, gamma_B,
            alpha_C, beta_C, gamma_C,
            alpha_D, beta_D, gamma_D,
            alpha_E, beta_E, gamma_E):
    """
    Equations 5.20 - 5.22 in KPark thesis (p. 71)
    """
    A = (alpha_A * theta_e**2 + beta_A * theta_e + gamma_A) * e_phi**4
    B = (alpha_B * theta_e**2 + beta_B * theta_e + gamma_B) * e_phi**3
    C = (alpha_C * theta_e**2 + beta_C * theta_e + gamma_C) * e_phi**2
    D = (alpha_D * theta_e**2 + beta_D * theta_e + gamma_D) * e_phi
    E = (alpha_E * theta_e**2 + beta_E * theta_e + gamma_E)

    return A + B + C + D + E


def Dtheta(e_phi, theta_e, A, B, C, D):
    """
    Mom Corrections for e6 (CLAS-NOTE 2003-005)
    """
    # e_phi = np.rad2deg(e_phi)
    first = (A+B*e_phi)*(np.cos(theta_e)/np.cos(e_phi))
    second = (C+D*e_phi)*np.sin(theta_e)
    return first + second


def Dpp(e_phi, theta_e, p, Bt, E, F, G, H):
    first = (E+F*e_phi)*(np.cos(theta_e)/np.cos(e_phi))
    second = (G+H*e_phi)*np.sin(theta_e)

    return (first + second)*(p/Bt)


@timeit
def read_file(file_name: str) -> pd.DataFrame:
    dtype = {
        "e_p": "float32",
        "e_theta": "float32",
        "e_phi": "float32",
        "p_p": "float32",
        "p_theta": "float32",
        "p_phi": "float32",
        "W_uncorr": "float32",
        "Q2_uncorr": "float32",
        "sector": "int8",
        "e_phi_center": "float32",
        "e_theta_calc": "float32",
        "delta_theta":   "float32",
        "p_p_calc":      "float32"
    }

    pyTable = csv.read_csv(
        file_name,
        read_options=csv.ReadOptions(use_threads=True),
        convert_options=csv.ConvertOptions(column_types=dtype)
    )
    df = pyTable.to_pandas(strings_to_categorical=True)

    df = df[(df.W_uncorr > 0.5) & (df.W_uncorr < 1.5)]
    df['e_theta'] = np.deg2rad(df.e_theta)
    df['e_phi_center'] = center_phi(df.e_phi, df.sector)
    df['e_phi'] = np.deg2rad(df.e_phi, dtype=np.float32)
    df['p_theta'] = np.deg2rad(df.p_theta, dtype=np.float32)
    df['p_phi'] = np.deg2rad(df.p_phi, dtype=np.float32)
    df['e_theta_calc'] = Theta_e_calc(df.p_theta)
    df['delta_theta'] = df['e_theta_calc']-df['e_theta']
    # df['w_corr'] = calc_W(df.e_p, df.e_theta, df.e_phi)
    # df['q2_corr'] = q2_calc(df.e_p, df.e_theta, df.e_phi)
    df['p_p_calc'] = Theta_e_calc(df.e_theta_calc)
    df.dropna(inplace=True)

    dtype_pd = {
        "e_p": "float16",
        "e_theta": "float16",
        "e_phi": "float16",
        "p_p": "float16",
        "p_theta": "float16",
        "p_phi": "float16",
        "W_uncorr": "float16",
        "Q2_uncorr": "float16",
        "sector": "int8",
        "e_phi_center": "float16",
        "e_theta_calc": "float16",
        "delta_theta":   "float16",
        "p_p_calc":      "float16"
    }
    df = df.astype(dtype_pd)

    return df


def get_min_max(x, past=0.1):
    a, b = np.min(x), np.max(x)
    if a <= 0:
        a, b = b, a

    return ((1-past)*a, (1+past)*b)


@timeit
def dtheta_vs_phi(sector_data: pd.DataFrame) -> Dict:
    outputs = dict()
    sec = int(np.mean(sector_data.sector))

    phi_step = 0.5
    phi_steps = np.arange(-30, 30, phi_step)
    theta_step = 1
    theta_steps = np.arange(13, 25, 1)

    for theta in theta_steps:
        phi_theta = []
        fig, ax = plt.subplots(figsize=(12, 9))
        dt_theta = sector_data[(np.rad2deg(sector_data.e_theta) > theta) & (
            (np.rad2deg(sector_data.e_theta)) <= theta+theta_step)]
        ax.hist2d(dt_theta.e_phi_center.to_numpy(),
                  dt_theta.delta_theta.to_numpy(), bins=250, range=[[-30, 30], [-0.01, 0.01]])

        for phi in phi_steps:
            dt = dt_theta[(dt_theta.e_phi_center > phi) & (
                dt_theta.e_phi_center <= phi+phi_step)]

            if len(dt) < 100:
                continue

            phi = np.mean(dt.e_phi_center.to_numpy())
            delta_t = np.mean(dt.delta_theta.to_numpy())

            try:
                fig2, ax2 = plt.subplots(figsize=(12, 9))
                yy, xx, _ = ax2.hist(dt.delta_theta.to_numpy(),
                                     bins=1000, range=(-0.05, 0.05))
                xx = (xx[1:]+xx[:-1])/2.0
                g_mod = BreitWignerModel()
                g_pars = g_mod.guess(yy, x=xx)
                g_out = g_mod.fit(yy, g_pars, x=xx)
                np.mean(dt.delta_theta.to_numpy())-g_out.best_values['center']
                xxs = np.linspace(-0.05, 0.05, 2000)
                ax2.plot(xxs, g_mod.eval(params=g_out.params, x=xxs))
                delta_t = g_out.best_values['center']
                yerr = g_out.best_values['sigma']
                # fig2.savefig(f"plots/slices/slice_{theta}_{phi}_{sec}.png")
                ax.errorbar(phi, delta_t,
                            yerr=yerr, fmt='.', c='red')
                del fig2
                del ax2
            except ValueError:
                pass

            phi_theta.append([phi, delta_t])

        x = np.transpose(phi_theta)[0]
        y = np.transpose(phi_theta)[1]

        n = 4
        z = np.polyfit(x, y, n)

        func = np.poly1d(z)
        xs = np.linspace(-30, 30, 2000)
        ax.plot(xs, func(xs), label=f'{z}')

        # popt, pcov = curve_fit(Dtheta, x, y, maxfev=3400)
        # ax.plot(xs, Dtheta(xs, *popt), label=f'{popt}')

        # popt, pcov = curve_fit(FitFunc, x, y, maxfev=3400)
        # ax.plot(xs, FitFunc(xs, *popt), label=f'{popt}')
        ax.set_ylim(-0.01, 0.01)
        ax.set_xlabel(f"$\phi$")
        ax.set_ylabel(f"$\Delta \\theta$")
        ax.legend()
        fig.savefig(f'plots/fit_sec_{sec}_theta_{theta}.png')
        outputs[f'sec_{sec}_theta_{theta}'] = tuple(z)
        del fig, ax

    return outputs


def slices(df: pd.DataFrame, num_points: int = 25) -> pd.DataFrame:
    grr = []
    for sec in range(1, 7):
        for deg in range(13, 20):
            df2 = df[(df.sector == sec) & (np.rad2deg(df.e_theta) >= deg)
                     & (np.rad2deg(df.e_theta) < deg+1)]
            if len(df2) < 2:
                pass
            phis = np.linspace(np.min(df2.e_phi_center),
                               np.max(df2.e_phi_center), num_points+1)
            for phi in range(0, num_points):
                phi_min = phis[phi]
                phi_max = phis[phi+1]
                data = df2[(df2.e_phi_center > phi_min) &
                           (df2.e_phi_center <= phi_max)]
                y, x = np.histogram(data['delta_theta'],
                                    bins=100, range=(-0.006, 0.006))
                y = y/np.max(y)
                x = (x[:-1]+x[1:])/2.0
                # plt.errorbar(x,y,fmt='o-',alpha=0.4)
                if not np.any(np.isnan(y)) and len(y) >= 50:
                    try:
                        mod = GaussianModel()
                        pars = mod.guess(y, x=x)
                        out = mod.fit(y, pars, x=x)

                        if np.abs(out.params['center']) < 0.2:
                            plt.plot(x, out.best_fit /
                                     np.max(out.best_fit), alpha=0.1)
                            grr.append(
                                [((phi_min+phi_max) / 2.0), out.params['center'], out.params['fwhm'], sec, deg])
                    except RuntimeError:
                        pass

    plt.show()
    grr = np.array(grr)
    return pd.DataFrame(data=grr, columns=["phi", "dtheta", "sigma", "sec", "min_deg"])


if __name__ == "__main__":
    print(np.__version__)
    if len(sys.argv) <= 1:
        exit()

    print(sys.argv[1])
    df = read_file(sys.argv[1])

    data_to_fit = [df[df.sector == sec] for sec in range(1, 2)]
    del df
    fit_data = dict()
    for sec, data in enumerate(data_to_fit):
        fit_data[sec] = dtheta_vs_phi(data)

    print(fit_data)
    # with Pool(6) as p:
    #     data_to_fit = p.map(dtheta_vs_phi, data_to_fit)

    # for sec in range(1, 7):
    #     fig, ax = plt.subplots(figsize=(12, 9))
    #     data = df[(df.sector == sec)]
    #     step_size = 1
    #     for phi in np.arange(min(np.degrees(data.e_phi)), max(np.degrees(data.e_phi)), step_size):
    #         dt = data[(np.degrees(data.e_phi) > phi) & (
    #             np.degrees(data.e_phi) <= phi+step_size)]
    #         if len(dt) < 10:
    #             continue
    #         ax.hist(dt.W_uncorr, bins=100, range=(0.8, 1.05), alpha=0.2)
    #     plt.savefig(f'plots/w_{sec}.png')
