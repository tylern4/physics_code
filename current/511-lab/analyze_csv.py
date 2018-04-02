#!/usr/bin/env python
from ROOT import TH2D, TH1D, TVector3, TLorentzVector, TFile
import numpy as np
import pandas as pd

#Setup beam 4 vector
BEAM = 4.81726  # Beam energy in GeV
MASS_P = 0.93827203
MASS_E = 0.000511
e_mu = TLorentzVector(0.0, 0.0, np.sqrt(BEAM * BEAM - MASS_E * MASS_E), BEAM)
p_target = TLorentzVector(0, 0, 0, MASS_P)
e_mu_prime_3 = TVector3()
e_mu_prime = TLorentzVector()

# A basic analysis program looking at
# a electron beam on a proton target

# Follow the TODO portions to get the analysis working properly


# TODO: Write functions to calculate W and Q^2
# Calcuating Q^2
# q^mu^2 = (e^mu - e^mu')^2 = -Q^2
def Q2(event):
    e_mu_prime_3.SetXYZ(event.px, event.py, event.pz)
    e_mu_prime.SetVectM(e_mu_prime_3, MASS_E)
    q_mu = (e_mu - e_mu_prime)
    return -q_mu.Mag2()


# Calcualting W
# gamma = e^mu - e^mu'
# P is proton target at rest
# Gotten from s channel [(gamma - P)^2 == s == w^2]
def W(event):
    e_mu_prime_3.SetXYZ(event.px, event.py, event.pz)
    e_mu_prime.SetVectM(e_mu_prime_3, MASS_E)
    q_mu = (e_mu - e_mu_prime)
    return (p_target + q_mu).Mag()


def analyze():
    fin = "/Users/tylern/physics_code/current/bin/511_lab_small.csv.gz"
    fout = "output.root"

    # TODO: Create the histograms you want to fill

    # Load chain from branch lab
    OutputFile = TFile(fout, "RECREATE")

    #Create 4 vectors for the scattered electron
    e_mu_prime_3 = TVector3()
    e_mu_prime = TLorentzVector()

    data = pd.read_csv(fin)

    data['px'] = data['p'] * data['cx']
    data['py'] = data['p'] * data['cy']
    data['pz'] = data['p'] * data['cz']

    data['Q2'] = data.apply(Q2, axis=1)
    data['W'] = data.apply(W, axis=1)
    print(data.head())
    #for event in data.itertuples():
    # We setup the scattered electron by first setting the momentum 3 vector
    #    e_mu_prime_3.SetXYZ(event.px, event.py, event.pz)
    #  And then adding a mass to the 3 vector
    #  ROOT will calculate the proper energy for
    #  the 4 vector for us this way.
    #  Check out (https://root.cern.ch/doc/master/classTLorentzVector.html)
    #  for reference
    #    e_mu_prime.SetVectM(e_mu_prime_3, MASS_E)
    #print(e_mu_prime.E())

    #TODO: Calculate W and Q^2 using the functions you created above and fill the histograms.

    #end stuff
    OutputFile.cd()
    #TODO: Write the histograms into the file
    OutputFile.Close()


if __name__ == "__main__":
    analyze()
