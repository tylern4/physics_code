import sys
import os
import glob

import ROOT
from ROOT import TLorentzVector, TVector3
from ROOT import gBenchmark, gROOT, std
gROOT.SetBatch(True)

from constants import *
from physics import *
import plots

import numpy as np
import pandas as pd

from pathos.multiprocessing import Pool, cpu_count, ProcessingPool
import multiprocessing as mp
try:
    import dill as pl
except:
    sys.exit()

import cppyy
cppyy.load_reflection_info("H10.so")


class datahandeler(object):
    """Datahandeler class"""

    def __init__(self, args):
        self.args = args
        self.args.ncore = (cpu_count(), self.args.ncore)[self.args.ncore > 0]
        if self.args.clean:
            self.clean()
        print("Starting datahandeler with %d cores" % self.args.ncore)

    def clean(self):
        print("Cleaning outputs in %s" % self.args.output)
        files = glob.glob(self.args.output + "*.pkl")
        for f in files:
            if os.path.exists(f):
                os.remove(f)

    def split_list(self):
        wanted_parts = self.args.ncore
        alist = glob.glob(self.args.input + '*.root')
        length = len(alist)
        return [alist[i * length // wanted_parts: (i + 1) * length // wanted_parts]
                for i in range(wanted_parts)]

    def _run(self, files):
        gBenchmark.Start("loop " + str(mp.current_process().pid))
        chain = ROOT.TChain('h10')
        h10 = cppyy.gbl.H10()
        for _f in files:
            chain.Add(_f)
        h10.loop(chain)
        W = np.array([_W for _W in h10.W_vec])
        Q2 = np.array([_Q2 for _Q2 in h10.Q2_vec])
        pl.dump(W, open(self.args.output + 'W_' +
                        str(mp.current_process().pid) + '.pkl', 'wb'), 2)
        pl.dump(Q2, open(self.args.output + 'Q2_' +
                         str(mp.current_process().pid) + '.pkl', 'wb'), 2)
        gBenchmark.Show("loop " + str(mp.current_process().pid))

    def run_mp(self):
        files = self.split_list()
        pool = ProcessingPool()
        out = pool.map(self._run, (files))
        pool.close()
        pool.join()

    def plot(self):
        plot = plots.plotting()
        plot.WvsQ2(self.W_Q2, output=self.args.output)
