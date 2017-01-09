#!/usr/bin/env python
import ROOT
import sys
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import pandas as pd
import argparse
from main import *
from constants import *
from physics import *


parser = argparse.ArgumentParser(description='Test Root Funtionality')
parser.add_argument('input', type=str, nargs='?')
parser.add_argument('output', type=str, nargs='?')

args = parser.parse_args()

if args.input is None:
	args.input = '.'

if args.output is None:
	args.output = '.'

directory = args.input + '/*.root'
chain_h10 = ROOT.TChain('h10')
chain_h10.Add(directory)

Q2 = np.array([])
W = np.array([])

for _e in chain_h10:
	if _e.id[0] is ID['ELECTRON']:
		e_mu_p = fourvec(_e.p[0],_e.cx[0],_e.cy[0],_e.cz[0],mass['ELECTRON'])
		Q2 = append(Q2,Q2_calc(e_mu,e_mu_p))
		W = append(W,W_calc(e_mu,e_mu_p))

fig = plt.figure(num=None, figsize=fig_size, dpi=200, facecolor='w', edgecolor='k')
plt.hist2d(W,Q2,bins=500,range=[[0,4],[0,10]],cmap=color_map)
plt.ylabel(r'$Q^{2}$ $(GeV^{2})$', fontsize=18)
plt.xlabel(r'$W (GeV)$', fontsize=20)
plt.colorbar()
plt.savefig(args.output + '/' + 'WvsQ2.pdf')

		


