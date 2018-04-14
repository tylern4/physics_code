#!/usr/bin/env python
from ROOT import TVector3, TLorentzVector
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

BEAM = 4.81726  # Beam energy in GeV
MASS_P = 0.93827203  # Mass in GeV
MASS_E = 0.000511  # Mass in GeV

Square = lambda x: np.square(x)
append = lambda _arr, _val: np.append(_arr, _val)

_e_mu = TLorentzVector(0.0, 0.0, np.sqrt(Square(BEAM) - Square(MASS_E)), BEAM)

_p_target = TLorentzVector(0, 0, 0, MASS_P)


def Q2_calc(e_mu_prime):
    """Retruns Q^2 value: q^mu^2 = (e^mu - e^mu')^2 = -Q^2."""
    _q_mu = (_e_mu - e_mu_prime)
    return -_q_mu.Mag2()


def W_calc(e_mu_prime):
    """Returns W: Gotten from s channel [(gamma - P)^2 == s == w^2], Sqrt[M_p^2 - Q^2 + 2 M_p gamma]."""
    _q_mu = (_e_mu - e_mu_prime)
    return (_p_target + _q_mu).Mag()


def analyze():
    fin = "511_lab_E_data.csv.gz"

    #Create 4 vectors for the scattered electron
    e_mu_prime_3 = TVector3()
    e_mu_prime = TLorentzVector()

    data = pd.read_csv(fin, sep=',')

    data['px'] = data['p'] * data['cx']
    data['py'] = data['p'] * data['cy']
    data['pz'] = data['p'] * data['cz']

    W = []
    Q2 = []

    for event in data.itertuples():
        e_mu_prime_3.SetXYZ(event.px, event.py, event.pz)
        e_mu_prime.SetVectM(e_mu_prime_3, MASS_E)
        W.append(W_calc(e_mu_prime))
        Q2.append(Q2_calc(e_mu_prime))

    # W
    output_file = "W.pdf"
    # Make the figure for plotting
    fig = plt.figure(
        num=None, figsize=(16, 9), dpi=200, facecolor='w', edgecolor='k')

    # Fill the histogram
    plt.hist(W, 500, normed=1, histtype='stepfilled', alpha=0.75, range=[0, 3])
    # Add labels
    plt.title("W")
    plt.xlabel("W (GeV)")
    plt.ylabel("Count")
    # Save to output_file name
    plt.savefig(output_file)

    # Q2
    output_file = "Q2.pdf"
    # Make the figure for plotting
    fig = plt.figure(
        num=None, figsize=(16, 9), dpi=200, facecolor='w', edgecolor='k')
    # Fill the histogram
    plt.hist(
        Q2, 500, normed=1, histtype='stepfilled', alpha=0.75, range=[0, 4])
    # Add labels
    plt.title("$Q^2$")
    plt.xlabel("$Q^2$ ($GeV^2$)")
    plt.ylabel("Count")
    # Save to output_file name
    plt.savefig(output_file)

    # WvsQ2
    output_file = "WvsQ2.pdf"
    # Make the figure for plotting
    fig = plt.figure(
        num=None, figsize=(16, 9), dpi=200, facecolor='w', edgecolor='k')
    # Fill the histogram
    plt.hist2d(W, Q2, bins=500, normed=True, range=[[0, 3], [0, 4]])
    # Add labels
    plt.title("W vs $Q^2$")
    plt.xlabel("W (GeV)")
    plt.ylabel("$Q^2$ ($Gev^2$)")
    # Save to output_file name
    plt.savefig(output_file)


if __name__ == "__main__":
    analyze()
