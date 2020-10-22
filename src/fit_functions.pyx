# distutils: language = c++
#cython: boundscheck=False
#cython: nonecheck=False
## cython: wraparound=False
#cython: infertypes=True
#cython: initializedcheck=False
#cython: cdivision=True
#cython: embedsignature=True
#cython: profile=True
#cython: binding=True


import numpy as np
cimport numpy as np
from scipy.special import erfc


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
