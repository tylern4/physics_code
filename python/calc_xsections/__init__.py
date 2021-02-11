from maid_interface import maid_2007_Npi as maid
import numpy as np
import pandas as pd
from pyarrow import csv
import time
from lmfit import Model
from lmfit.models import *
from scipy.special import erfc
from scipy import stats
import boost_histogram as bh

ENERGY = 4.81726

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


def luminosity():
    # Q_tot = 15623.89E-6  # mC E-6 -> C
    Q_tot = 31426.514045353397E-7  # 02/07/2021
    # Q_tot = 2906.9858369E-6 # 2nd attempt
    # Q_tot = 2742.5087E-6 # 1st attempt
    l = 5  # cm
    rho = 0.0708  # g/cm3
    Avigadro = 6.022E23  # mol^−1
    qe = 1.602E-19  # C
    MH = 1.007  # g/mol
    conv_cm2_to_fm2 = 1E-39  # From wolfram alpha
    conv_cm2_to_mubarn = 1E-30  # From wolfram alpha

    # * conv_cm2_to_fm2
    return (Q_tot*l*rho*Avigadro)/(qe*MH) * conv_cm2_to_mubarn


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


def hist_data(data, density=True, bins=10):
    data_y, data_x = bh.numpy.histogram(
        data.phi.to_numpy(), bins=bins, range=(0, 2 * np.pi), density=density, threads=4)
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
    F = mc_rec_y/thrown_y
    error = (thrown_y-mc_rec_y)*mc_rec_y
    error = error / (thrown_y**3)
    error = np.sqrt(error)
    error = (error/F)
    # error_bar = np.sqrt(y*error**2 + stats.sem(y)**2)
    # error_bar = np.sqrt(y*error**2 + np.std(y)**2)
    error_bar = np.sqrt(y*error**2 + np.std(y)**2)

    return error_bar


def plot_maid_model(ax, w, q2, theta, xs):
    # Get bin centers
    _w = (w.left + w.right) / 2.0
    _q2 = (q2.left + q2.right) / 2.0
    _theta = (theta.left + theta.right) / 2.0
    # Get the cross section values from maid
    crossSections = get_maid_values(xs, _w, _q2, _theta)
    # _ax = ax.twinx()
    ax.plot(xs, crossSections, c='r', linestyle='dotted')
    # ax.set_ylim(bottom=0, top=np.max(crossSections)*1.5)

    return ax


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
            linewidth=2.0, c=color, label=f'{name}', alpha=0.2)

    # Get uncertinty and plot between 2 sigmas
    # dely = out.eval_uncertainty(sigma=3, x=xs)
    # ax.fill_between(xs,
    #                 out.eval(params=out.params, x=xs)-dely,
    #                 out.eval(params=out.params, x=xs)+dely,
    #                 color=color, alpha=0.1,
    #                 label='$3 \sigma$')

    return out
