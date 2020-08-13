#!/usr/bin/env python
from ROOT import TVector3, TLorentzVector
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

BEAM = 4.81726  # Beam energy in GeV
MASS_P = 0.93827203  # Mass in GeV
MASS_E = 0.000511  # Mass in GeV
MASS_PIP = 0.13957018

Square = lambda x: np.square(x)
append = lambda _arr, _val: np.append(_arr, _val)

_e_mu = TLorentzVector(0.0, 0.0, np.sqrt(Square(BEAM) - Square(MASS_E)), BEAM)

_p_target = TLorentzVector(0, 0, 0, MASS_P)


def MM_calc(e_mu_prime, pip_mu_prime):
    _q_mu = (_e_mu - e_mu_prime)
    return (_q_mu + _p_target - pip_mu_prime).Mag()


def analyze():
    fin = "511_lab_E_PIP_data.csv"

    #Create 4 vectors for the scattered electron
    e_mu_prime_3 = TVector3()
    e_mu_prime = TLorentzVector()

    pip_mu_prime_3 = TVector3()
    pip_mu_prime = TLorentzVector()

    data = pd.read_csv(fin, sep=',')

    data['e_px'] = data['e_p'] * data['e_cx']
    data['e_py'] = data['e_p'] * data['e_cy']
    data['e_pz'] = data['e_p'] * data['e_cz']
    data['pip_px'] = data['pip_p'] * data['pip_cx']
    data['pip_py'] = data['pip_p'] * data['pip_cy']
    data['pip_pz'] = data['pip_p'] * data['pip_cz']

    mm = []

    for event in data.itertuples():
        e_mu_prime_3.SetXYZ(event.e_px, event.e_py, event.e_pz)
        e_mu_prime.SetVectM(e_mu_prime_3, MASS_E)
        pip_mu_prime_3.SetXYZ(event.pip_px, event.pip_py, event.pip_pz)
        pip_mu_prime.SetVectM(pip_mu_prime_3, MASS_PIP)
        mm.append(MM_calc(e_mu_prime, pip_mu_prime))

    # W
    output_file = "mm.pdf"
    # Make the figure for plotting
    fig = plt.figure(
        num=None, figsize=(16, 9), dpi=200, facecolor='w', edgecolor='k')

    # Fill the histogram
    plt.hist(W, 500, normed=1, histtype='stepfilled', alpha=0.75, range=[0, 3])
    # Add labels
    plt.title("Missing Mass")
    plt.xlabel("Missing Mass (GeV)")
    plt.ylabel("Count")
    # Save to output_file name
    plt.savefig(output_file)


if __name__ == "__main__":
    analyze()
