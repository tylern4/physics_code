#!/usr/bin/env python
# coding: utf-8

# In[1]:


import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')

MASS_P = 0.93827203
MASS_E = 0.000511
E0 = 4.81726


# In[2]:


df = pd.read_csv("/Users/tylern/Desktop/show/mom_corr.dat")


# In[3]:


df["theta_e_measured"] = np.arccos(
    df.e_pz/np.sqrt(df.e_px**2 + df.e_py**2 + df.e_pz**2))
df["theta_p_measured"] = np.arccos(
    df.p_pz/np.sqrt(df.p_px**2 + df.p_py**2 + df.p_pz**2))
df["phi_e"] = np.arctan2(df.e_px/np.sqrt(df.e_px**2 + df.e_py**2 + df.e_pz**2),
                         df.e_py/np.sqrt(df.e_px**2 + df.e_py**2 + df.e_pz**2))
df["theta_e_calc"] = 2 * \
    np.arctan(MASS_P / ((E0 + MASS_P) * np.tan(df.theta_p_measured)))
df["delta_theta"] = df.theta_e_calc - df.theta_p_measured


# In[4]:


df = df[df["Q2_uncorr"] < 1.4]
df = df[df["W_uncorr"] <= 1.05]
df = df[df["W_uncorr"] >= 0.7]
df = df[np.degrees(df["theta_p_measured"]) > 35]


# In[5]:


df.head()


# In[6]:


for i in range(1, 7):
    data = df[(df.sector == i)]
    step_size = 1
    for phi in np.arange(min(np.degrees(data.phi_e)), max(np.degrees(data.phi_e)), step_size):
        dt = data[(np.degrees(data.phi_e) > phi) & (
            np.degrees(data.phi_e) <= phi+step_size)]
        if len(dt) < 10:
            continue
        plt.hist(dt.W_uncorr, bins=100, range=(0.8, 1.05))
        # plt.show()
        #mean = np.mean(dt.delta_theta)
        #std = np.std(dt.delta_theta)
        #plt.hist(dt.delta_theta, bins=100, range=(mean-3*std, mean+3*std))
        # plt.show()
    plt.show()


# In[58]:


for i in range(1, 7):
    data = df[(df.sector == i)]
    step_size = 1
    phi_theta = []
    for phi in np.arange(min(np.degrees(data.phi_e)), max(np.degrees(data.phi_e)), step_size):
        dt = data[(np.degrees(data.phi_e) > phi) & (
            np.degrees(data.phi_e) <= phi+step_size)]
        if len(dt) < 10:
            continue
        plt.errorbar(np.mean(dt["phi_e"]), np.mean(
            dt["delta_theta"]), yerr=np.std(dt["delta_theta"]), fmt='o', c='blue')
        phi_theta.append(np.array([np.mean(dt["phi_e"]), np.mean(
            dt["delta_theta"]), np.std(dt["delta_theta"])]))
    x = np.transpose(phi_theta)[0]
    y = np.transpose(phi_theta)[1]
    yerr = np.transpose(phi_theta)[2]
    n = 4
    z = np.polyfit(x, y, n, w=1/(2*yerr))
    func = np.poly1d(z)
    plt.plot(x, func(x), label=f'{z}')
    plt.legend()
    plt.show()


# In[ ]:
