#!/usr/bin/env python
import numpy as np
from numpy.random import normal
import matplotlib.pyplot as plt
import scipy.stats as stats


def func(num):
    a = 0.5
    b = 2.0
    c = 1.0
    x = np.random.uniform(-1.0, 0.6, num)
    y = a*x**2 + b*x + c
    return y


num = 100000
data1 = normal(loc=0.98, scale=0.1, size=200000)
data2 = func(700000)
#data2 = normal(loc=1.02, scale=0.3, size=150000)


plt.hist(data1, bins=100)
plt.hist(data2, bins=100)

plt.show()

data = np.concatenate([data1, data2])
plt.hist(data, bins=200, density=True)

# lets try the normal distribution first
m, s = stats.norm.fit(data)

m = data.mean()
s = data.std()
# get mean and standard deviation
# now get theoretical values in our interval
lnspc = np.linspace(0.5, 1.2, 500)
pdf_g = stats.norm.pdf(lnspc, m, s)
plt.plot(lnspc, pdf_g, label="Norm")  # plot it

plt.show()

for d in data:
    print(d)
