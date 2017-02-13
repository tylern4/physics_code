import ROOT
from ROOT import gROOT, gBenchmark
from physics import *
import numpy as np

import cppyy
cppyy.load_reflection_info("Physics.so")
physics = cppyy.gbl.Physics()

gROOT.SetBatch(True)
chain = ROOT.TChain('h10')
num_files = chain.Add('/Users/tylern/Desktop/dumb/r22855_00.root')

Q2 = np.array([])
W = np.array([])

# looping over events
gBenchmark.Start('Python calcs')
for evnt in chain:
    if evnt.id[0] is get_id['ELECTRON']:
        e_mu_p = fourvec(evnt.p[0], evnt.cx[0],
                         evnt.cy[0], evnt.cz[0], get_mass('ELECTRON'))
        Q2 = append(Q2, Q2_calc(e_mu, e_mu_p))
        W = append(W, W_calc(e_mu, e_mu_p))
gBenchmark.Show('Python calcs')

# looping over events
gBenchmark.Start('CPP calcs')
for evnt in chain:
    if evnt.id[0] is get_id['ELECTRON']:
        e_mu_p = fourvec(evnt.p[0], evnt.cx[0],
                         evnt.cy[0], evnt.cz[0], get_mass('ELECTRON'))
        Q2 = append(Q2, physics.Q2_calc(e_mu, e_mu_p))
        W = append(W, physics.W_calc(e_mu, e_mu_p))
gBenchmark.Show('CPP calcs')
