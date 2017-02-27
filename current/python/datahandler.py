import sys
import os
import glob

import ROOT
from ROOT import TLorentzVector, TVector3, TH2D, TFile
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

import sys


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
        hist_obj = h10.loop_test(chain)

        # h10 contains values from the c++ class in a ROOT vector form which
        # can be iterated over
        # W_Q2 = np.array([_W for _W in h10.W_vec], [_Q2 for _Q2 in h10.Q2_vec])
        # Q2 = np.array()

        # TODO: Take vector values and send them to be plotted
        # At this point I'm dumping the numpy vectors to pickles to be able to
        # graph them later.
        # pl.dump(W, open(self.args.output + 'W_' +
        #                str(mp.current_process().pid) + '.pkl', 'wb'), 2)
        # pl.dump(Q2, open(self.args.output + 'Q2_' +
        #                 str(mp.current_process().pid) + '.pkl', 'wb'), 2)
        gBenchmark.Show("loop " + str(mp.current_process().pid))

        histograms = {'WvsQ2_hist': hist_obj.WvsQ2_hist}
        return histograms

    def run_map(self):
        """Maps function to run on multiple cores"""
        # Split input into chunks for processing
        files = self.split_list()
        # Make processing pool
        pool = ProcessingPool()
        # Map processing to _run function
        self.output = pool.map(self._run, (files))
        # Close and join pool
        pool.close()
        pool.join()

    def run_reduce(self):
        """Reduce function to add all histograms at the end and save"""
        file = TFile(self.args.output + "test.root", "RECREATE")
        WvsQ2_hist = TH2D("WvsQ2_hist", "W vs Q^{2}", 500, 0, 3.25,
                          500, 0, 10)

        for _h in self.output:
            WvsQ2_hist.Add(_h['WvsQ2_hist'])

        WvsQ2_hist.Write()
        file.Write()
        file.Close()

    def plot(self):
        plot = plots.plotting()
        plot.WvsQ2(self.W_Q2, output=self.args.output)
