#!/usr/local/bin/ipython
import numpy as np
import pandas as pd
import xray
import scipy
from netCDF4 import num2date
import matplotlib.pyplot as plt
from tqdm import *
from scipy.optimize import curve_fit

np.set_printoptions(precision=4, suppress=True)
plt.rc('text', usetex=True)
num_bins = 500
fig_size = (16, 9)

lines = [line.rstrip('\n') for line in open('/home/tylern/physics_code/python/v2_all.lis')]
bad = []
for line in tqdm(lines):
    try:
        line = line.replace('root','h5')
        store = pd.HDFStore(line)
    except:
        bad.append(line.replace('h5','root'))
print(bad)
df = store['df']
xr = xray.Dataset.from_dataframe(df)

def pi_momentum_B(dataset):
    pi_p , pi_b = [],[]
    ids = dataset['id'].to_series()
    P = dataset['p'].to_series()
    B = dataset['b'].to_series()
    
    for i in tqdm(range(len(ids))):
        for j in range(len(ids[i])):
            #if ids[i][j] != 0 and ids[i][j] != 11:
            if ids[i][j] == 211:
                pi_p.append(P[i][j])
                pi_b.append(B[i][j])
    return np.array(pi_p),np.array(pi_b)



fig = plt.figure(num=None, figsize=fig_size, dpi=200, facecolor='w', edgecolor='k')
plt.hist2d(pip_p, pip_b, bins=num_bins,cmap='viridis',range=[[0,1.2],[0.1,1.2]])
plt.colorbar()
plt.show()

def charge_momentum_B(dataset):
    pos_p , pos_b = [],[]
    neg_p , neg_b = [],[]
    charge = dataset['q'].to_series()
    P = dataset['p'].to_series()
    B = dataset['b'].to_series()
    ids = dataset['id'].to_series()
    for i in tqdm(range(len(charge))):
        for j in range(len(charge[i])):
            if charge[i][j] == 1:
                pos_p.append(P[i][j])
                pos_b.append(B[i][j])
            elif charge[i][j] == -1 and ids[i][j] != 11:
                neg_p.append(P[i][j])
                neg_b.append(B[i][j])
    return np.array(pos_p),np.array(pos_b),np.array(neg_p),np.array(neg_b)

p_pos,b_pos,p_neg,b_neg = charge_momentum_B(xr)
fig = plt.figure(num=None, figsize=fig_size, dpi=200, facecolor='w', edgecolor='k')
plt.hist2d(p_pos, b_pos, bins=num_bins,cmap='viridis',range=[[0,1.2],[0.1,1.2]])
plt.colorbar()
plt.show()
fig2 = plt.figure(num=None, figsize=fig_size, dpi=200, facecolor='w', edgecolor='k')
plt.hist2d(p_neg, b_neg, bins=num_bins,cmap='viridis',range=[[0,1.2],[0.1,1.2]])
plt.colorbar()

MM = xr['MM'].to_series()

def poly(x, c1, c2, c3):
    return c1*x*x + c2*x + c3

def gaussian(x, mu, sig, const):
    return const * 1/(sig*np.sqrt(2*np.pi)) * np.exp(-(x - mu)**2 / 2*sig**2)

def gaus_gaus(x, mu, sig, const, mu1, sig2, const3):
    return gaussian(x, mu1, sig2, const3) + gaussian(x, mu, sig, const)


hist, bin_edges = np.histogram(MM, bins=num_bins, range=(0,3))
xdata = 0.5*(bin_edges[1:]+bin_edges[:-1])
ydata = hist

#x0 = np.array([0.9447,66.3717,70871.596,1,1,1])

popt_1, pcov_1 = curve_fit(gaussian, xdata, ydata,maxfev =200000)
x0 = np.array([0.9447,50,1,popt_1[0],popt_1[1],popt_1[2]])
popt_1, pcov_1 = curve_fit(gaus_gaus, xdata, ydata, maxfev =200000,p0=x0)
perr_1 = np.sqrt(np.diag(pcov_1))


gaus_params = np.array([popt_1[3],popt_1[4],popt_1[5]])
fig1 = plt.figure(num=None, figsize=fig_size, dpi=200, facecolor='w', edgecolor='k')
plt.hist(MM, bins=num_bins, range=[0,3], histtype=u'stepfilled', facecolor='g', alpha=0.45)
plt.plot(xdata,gaus_gaus(xdata,*popt_1),'b--', lw=4)
plt.plot(xdata,gaussian(xdata,*popt_1[3:]),'r-.', lw=4)
plt.show()
#print popt_1
mass = popt_1[0]
sig = perr_1[0]
plt.show()