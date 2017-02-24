#!/usr/bin/env python
import ROOT
from ROOT import gROOT, gBenchmark
from pathos.multiprocessing import Pool, cpu_count, ProcessingPool
import numpy as np

from physics import *

# Load in the functions written by me in c++
import cppyy
# cppyy.load_reflection_info("Physics.so")
# physics = cppyy.gbl.Physics()
cppyy.load_reflection_info("H10.so")

# Sets root to batch mode to eliminate any pop up windows
gROOT.SetBatch(True)
# Make a chain to load in files, looking at the h10 branch in each file
chain = ROOT.TChain('h10')
# Load the files, Can also use glob to get all files in the folder.
num_files = chain.Add('~/Desktop/dumb/*.root')

# Make the Q^2 and W arrays that will be filled.
Q2 = np.array([])
W = np.array([])
evnt_list = []

import timeit


def calc(evnt):
    if evnt.id[0] is get_id['ELECTRON']:
        e_mu_p = fourvec(evnt.p[0], evnt.cx[0],
                         evnt.cy[0], evnt.cz[0], get_mass('ELECTRON'))
        Q2 = append(Q2, Q2_calc(e_mu, e_mu_p))
        W = append(W, W_calc(e_mu, e_mu_p))
        print(W_calc(e_mu, e_mu_p))


def load(chain):
    h10 = cppyy.gbl.H10()
    chain = h10.pass_chain(chain)
    return chain
    '''
    if evnt.id[0] is get_id['ELECTRON']:
        e_mu_p = fourvec(evnt.p[0], evnt.cx[0],
                         evnt.cy[0], evnt.cz[0], get_mass('ELECTRON'))
        Q2 = append(Q2, Q2_calc(e_mu, e_mu_p))
        W = append(W, W_calc(e_mu, e_mu_p))
        print(W_calc(e_mu, e_mu_p))
'''
pool = ProcessingPool()
c = pool.map(load, (chain))
pool.map(calc, (c))

print(W)
''''
# Start of ROOT based benchmark for testing
    gBenchmark.Start('Python calcs')
    # Loop over all events that have been added to the chain
    for evnt in chain:
        # If the events ID == ELECTRON then calcuate W and Q^2
        if evnt.id[0] is get_id['ELECTRON']:
            # Make electron four vector
            e_mu_p = fourvec(evnt.p[0], evnt.cx[0],
                             evnt.cy[0], evnt.cz[0], get_mass('ELECTRON'))
            # Calcuate W and Q^2 and append them to the arrays
            # Calcuations for W and Q^2 found in physics.py
            Q2 = append(Q2, Q2_calc(e_mu, e_mu_p))
            W = append(W, W_calc(e_mu, e_mu_p))
    gBenchmark.Show('Python calcs')
'''
'''
# Reset everying to do c++ based calculations to compare
# del W
# del Q2

# Make the Q^2 and W arrays that will be filled.
Q2 = np.array([])
W = np.array([])

# Start of ROOT based benchmark for testing
gBenchmark.Start('CPP calcs')
# Loop over all events that have been added to the chain
for evnt in chain:
    # If the events ID == ELECTRON then calcuate W and Q^2
    if evnt.id[0] is get_id['ELECTRON']:
        # Make electron four vector
        e_mu_p = fourvec(evnt.p[0], evnt.cx[0],
                         evnt.cy[0], evnt.cz[0], get_mass('ELECTRON'))
        # Calcuate W and Q^2 and append them to the arrays
        # Calcuations for W and Q^2 found in Physics.h and compiled with cppyy
        Q2 = append(Q2, physics.Q2_calc(e_mu, e_mu_p))
        W = append(W, physics.W_calc(e_mu, e_mu_p))
gBenchmark.Show('CPP calcs')
'''
