#!/usr/bin/env python
import numpy as np
from numpy.random import normal
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def func(x, a, b, c):
    y = a*x**2 + b*x + c
    y += normal(loc=1.0, scale=10.0, size=len(x))
    yerr = normal(loc=30.0, scale=10.0, size=len(x))
    return y, yerr


def fit_func(x, a, b, c):
    y = a*x**2 + b*x + c
    return y


xs = np.linspace(0, 50, 200)
x = np.random.randint(low=0, high=50, size=200)
a = -0.3
b = 0.2
c = 1000
y, yerr = func(x, a, b, c)
coef = np.polyfit(x, y, 2)
poly = np.poly1d(coef)

popt, pcov = curve_fit(fit_func, x, y)

plt.plot(xs, poly(xs), 'b.',
         label=f'a {coef[0]:0.2f} b {coef[1]:0.2f} c {coef[2]:0.2f}')
plt.plot(xs, fit_func(xs, *popt), 'r-',
         label=f'a {popt[0]:0.2f} b {popt[1]:0.2f} c {popt[2]:0.2f}')
plt.errorbar(x, y, yerr=yerr, fmt='.')
plt.legend()
# plt.show()


print("x,y,yerr")
for a, b, c in zip(x, y, yerr):
    print(a, b, c, sep=",")
