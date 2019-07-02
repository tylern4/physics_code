#!/usr/bin/env python
import ROOT
from ROOT import gROOT, gBenchmark
import numpy as np
import matplotlib
matplotlib.use('pgf')
import matplotlib.pyplot as plt
import argparse
import sys
from physics import *

parser = argparse.ArgumentParser(description="Python ROOT program")
parser.add_argument('input', type=str,
                    help="Input directory for *.root files")
parser.add_argument('output', type=str, nargs='?',
                    help="Output pdf file", default='./test.pdf')

# Print help if there aren't enough arguments
if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()

args = parser.parse_args()

if args.input[-1] != '/':
    args.input = args.input + '/'
if args.output[-4:] != '.pdf':
    args.output = args.output + '.pdf'

gROOT.SetBatch(True)
chain = ROOT.TChain('h10')
num_files = chain.Add(args.input + "*.root")

Q2 = []
W = []

for evnt in chain:
    if evnt.id[0] is get_id['ELECTRON']:
        e_mu_p = fourvec(evnt.p[0], evnt.cx[0], evnt.cy[0],
                         evnt.cz[0], get_mass('ELECTRON'))
        Q2.append(Q2_calc(e_mu, e_mu_p))
        W.append(W_calc(e_mu, e_mu_p))

fig = plt.figure(num=None, figsize=(16, 9), dpi=200,
                 facecolor='w', edgecolor='k')
plt.hist2d(W, Q2, bins=500, range=[[0, 3.14], [0, 10]])
plt.title("$W vs Q^{2}$")
plt.xlabel("$W (GeV)$")
plt.ylabel("$Q^{2} (GeV^{2})$")
plt.colorbar()
plt.savefig(args.output)
