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
        """Cleans the output directory of pickle files before saving new files"""
        print("Cleaning outputs in %s" % self.args.output)
        files = glob.glob(self.args.output + "*.pkl")
        for f in files:
            if os.path.exists(f):
                os.remove(f)

    def split_list(self):
        """Split the list of input files into equal chunks to be processed by the pool."""
        wanted_parts = self.args.ncore
        alist = glob.glob(self.args.input + '*.root')
        length = len(alist)
        return [alist[i * length // wanted_parts: (i + 1) * length // wanted_parts]
                for i in range(wanted_parts)]

    def _run(self, files):
        """Mapped run funtion"""
        # Start benchmark for each loop based on pid
        gBenchmark.Start("loop " + str(mp.current_process().pid))
        # Make a chain to load in files, looking at the h10 branch in each file
        chain = ROOT.TChain('h10')
        # Crete h10 object (Can be found in H10.h
        h10 = cppyy.gbl.H10()
        # For each file in the list add it to the chain to be processed
        for _f in files:
            chain.Add(_f)
        # Pass the chain to the h10 loop function (c++ based)
        h10.loop(chain)

        # h10 contains values from the c++ class in a ROOT vector form which
        # can be iterated over
        W = np.array([_W for _W in h10.W_vec])
        Q2 = np.array([_Q2 for _Q2 in h10.Q2_vec])
        # At this point I'm dumping the numpy vectors to pickles to be able to
        # graph them later.
        # My hope would be to return them to be plotted later.
        pl.dump(W, open(self.args.output + 'W_' +
                        str(mp.current_process().pid) + '.pkl', 'wb'), 2)
        pl.dump(Q2, open(self.args.output + 'Q2_' +
                         str(mp.current_process().pid) + '.pkl', 'wb'), 2)
        gBenchmark.Show("loop " + str(mp.current_process().pid))

    def run_mp(self):
        # Split input into
        files = self.split_list()
        pool = ProcessingPool()
        out = pool.map(self._run, (files))
        pool.close()
        pool.join()

    def plot(self):
        plot = plots.plotting()
        plot.WvsQ2(self.W_Q2, output=self.args.output)
