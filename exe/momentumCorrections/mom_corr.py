#!/usr/bin/env python

import pyarrow  # noqa
pyarrow.set_cpu_count(8)  # noqa
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


def dtheta_vs_phi(sector_data: pd.DataFrame, directory: str = ".") -> Dict:
    outputs = dict()
    sec = int(np.mean(sector_data.sector))

    phi_step = 0.5
    phi_steps = np.arange(-30, 30, phi_step)
    theta_step = 1
    theta_steps = np.arange(13, 45, theta_step)

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

        if len(phi_theta) > 3:
            x = np.transpose(phi_theta)[0]
            y = np.transpose(phi_theta)[1]

            n = 4
            z = np.polyfit(x, y, n)
            func = np.poly1d(z)
            xs = np.linspace(-30, 30, 2000)
            ax.plot(xs, func(xs), label=f'{z}')
            ax.legend()
            outputs[f'theta_{theta}'] = tuple(z)

        # popt, pcov = curve_fit(Dtheta, x, y, maxfev=3400)
        # ax.plot(xs, Dtheta(xs, *popt), label=f'{popt}')

        # popt, pcov = curve_fit(FitFunc, x, y, maxfev=3400)
        # ax.plot(xs, FitFunc(xs, *popt), label=f'{popt}')
        ax.set_ylim(-0.01, 0.01)
        ax.set_xlabel(f"$\phi$")
        ax.set_ylabel(f"$\Delta \\theta$")
        fig.savefig(f'{directory}/fit_sec_{sec}_theta_{theta}.png')

        del fig, ax

    return outputs


def second_fit(sector, fit_ABCDE, directory):
    theta_step = 1
    theta_steps = np.arange(13, 45, theta_step)
    thetas = []
    AS = []
    BS = []
    CS = []
    DS = []
    ES = []

    for theta in theta_steps:
        if f'theta_{theta}' in fit_ABCDE:
            thetas.append(theta)
            AS.append(fit_ABCDE[f'theta_{theta}'][0])
            BS.append(fit_ABCDE[f'theta_{theta}'][1])
            CS.append(fit_ABCDE[f'theta_{theta}'][2])
            DS.append(fit_ABCDE[f'theta_{theta}'][3])
            ES.append(fit_ABCDE[f'theta_{theta}'][4])

    thetas = np.array(thetas)
    AS = np.array(AS)
    BS = np.array(BS)
    CS = np.array(CS)
    DS = np.array(DS)
    ES = np.array(ES)

    xs = np.linspace(10, 30, 1000)
    n = 4
    fig, ax = plt.subplots(figsize=(12, 9))
    ax.set_ylim(-0.005, 0.005)

    ax.scatter(thetas, 10E5*AS, label="A", alpha=0.5, marker='o')
    z = np.polyfit(thetas, AS, n)
    func = np.poly1d(z)
    ax.plot(xs, 10E5*func(xs), label=f'{z}')

    ax.scatter(thetas, 10E5*BS, label="B", alpha=0.5, marker='s')
    z = np.polyfit(thetas, BS, n)
    func = np.poly1d(z)
    ax.plot(xs, 10E5*func(xs), label=f'{z}')

    ax.scatter(thetas, 10E5*CS, label="C", alpha=0.5, marker='^')
    z = np.polyfit(thetas, CS, n)
    func = np.poly1d(z)
    ax.plot(xs, 10E5*func(xs), label=f'{z}')

    ax.scatter(thetas, 10E5*DS, label="D", alpha=0.5, marker='v')
    z = np.polyfit(thetas, DS, n)
    func = np.poly1d(z)
    ax.plot(xs, 10E5*func(xs), label=f'{z}')

    ax.scatter(thetas, 10E5*ES, label="E", alpha=0.5, marker='>')
    z = np.polyfit(thetas, ES, n)
    func = np.poly1d(z)
    ax.plot(xs, 10E5*func(xs), label=f'{z}')

    ax.legend()
    fig.savefig(f"{directory}/fit_{sector}.png")


if __name__ == "__main__":
    print(np.__version__)
    if len(sys.argv) <= 1:
        exit()

    df = read_file(sys.argv[1])
    if len(sys.argv) == 3:
        output_folder = sys.argv[2]
    else:
        output_folder = "."

    data_to_fit = [df[df.sector == sec] for sec in range(1, 7)]
    output_folders = [output_folder for sec in range(1, 7)]
    del df

    ts = time.time()
    with Pool(6) as p:
        fit_ABCDE = p.starmap(
            dtheta_vs_phi, zip(data_to_fit, output_folders))

    te = time.time()
    print(f'First fits {clock_to_time(te - ts)}')

    ts = time.time()
    for sec, fit in enumerate(fit_ABCDE):
        second_fit(sec, fit, output_folder)

    te = time.time()
    print(f'Second fits {clock_to_time(te - ts)}')

    # with Pool(6) as p:
    #     fit_ABCD = p.starmap(
    #         dtheta_vs_phi, zip(data_to_fit, output_folders))
