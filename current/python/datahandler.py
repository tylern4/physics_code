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

# from pathos.multiprocessing import Pool, cpu_count, ProcessingPool
from multiprocessing import Pool, cpu_count
import multiprocessing as mp

import cppyy
cppyy.load_reflection_info("H10.so")

import sys


def _run(files):
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
    hist_obj = h10.dataHandeler(chain)

    gBenchmark.Show("loop " + str(mp.current_process().pid))
    histograms = {'WvsQ2_hist': hist_obj.WvsQ2_hist,
                  'W_hist': hist_obj.W_hist,
                  'Q2_hist': hist_obj.Q2_hist,
                  'E_prime_hist': hist_obj.E_prime_hist,
                  'Q2_vs_xb': hist_obj.Q2_vs_xb,
                  'WvsQ2_proton': hist_obj.WvsQ2_proton,
                  'W_proton': hist_obj.W_proton,
                  'Q2_proton': hist_obj.Q2_proton,
                  'WvsQ2_pion': hist_obj.WvsQ2_pion,
                  'W_pion': hist_obj.W_pion,
                  'Q2_pion': hist_obj.Q2_pion,
                  'WvsQ2_single_pi': hist_obj.WvsQ2_single_pi,
                  'W_single_pi': hist_obj.W_single_pi,
                  'Q2_single_pi': hist_obj.Q2_single_pi,
                  'WvsQ2_single_proton': hist_obj.WvsQ2_single_proton,
                  'W_single_proton': hist_obj.W_single_proton,
                  'Q2_single_proton': hist_obj.Q2_single_proton,
                  'MomVsBeta_hist': hist_obj.MomVsBeta_hist,
                  'MomVsBeta_hist_pos': hist_obj.MomVsBeta_hist_pos,
                  'MomVsBeta_hist_neg': hist_obj.MomVsBeta_hist_neg,
                  'Mom': hist_obj.Mom,
                  'Energy_hist': hist_obj.Energy_hist,
                  'MomVsBeta_proton_ID': hist_obj.MomVsBeta_proton_ID,
                  'MomVsBeta_Pi_ID': hist_obj.MomVsBeta_Pi_ID,
                  'MomVsBeta_proton_Pi_ID': hist_obj.MomVsBeta_proton_Pi_ID,
                  'Mass': hist_obj.Mass,
                  'Missing_Mass': hist_obj.Missing_Mass,
                  'Missing_Mass_square': hist_obj.Missing_Mass_square,
                  'delta_t_mass_P': hist_obj.delta_t_mass_P,
                  'delta_t_mass_P_PID': hist_obj.delta_t_mass_P_PID,
                  'delta_t_mass_PIP': hist_obj.delta_t_mass_PIP,
                  'delta_t_mass_PIP_PID': hist_obj.delta_t_mass_PIP_PID,
                  'EC_sampling_fraction': hist_obj.EC_sampling_fraction,
                  'fid_hist': hist_obj.fid_hist,
                  'fid_sec_hist': hist_obj.fid_sec_hist
                  }

    return histograms


def _run_skim(root_file_name):
    """Mapped skim funtion"""
    # Make a chain to load in files, looking at the h10 branch in each file
    chain = ROOT.TChain('h10')
    # Crete h10 object (Can be found in H10.h
    h10 = cppyy.gbl.H10()
    out_file_name = root_file_name.replace(
        "/root", "/skim").replace(".root", "_skim.root")
    h10.Skim(root_file_name, out_file_name)


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

    def run_map(self):
        """Maps function to run on multiple cores"""
        # Split input into chunks for processing
        files = self.split_list()
        # Make processing pool
        pool = Pool(processes=self.args.ncore)
        # Map processing to _run function
        self.output = pool.map(_run, files)
        # Close and join pool
        pool.close()
        pool.join()

    def run_reduce(self):
        """Reduce function to add all histograms at the end and save"""
        from histogram import add_and_save
        root_file = TFile(self.args.output, "RECREATE")
        add_and_save(self.output, root_file)

        root_file.Write()
        root_file.Close()

    def run_skim(self):
        """Maps function to skim on multiple cores"""
        # Split input into chunks for processin
        skim_files = glob.glob(self.args.input + "*.root")
        # Make processing pool
        pool = Pool(processes=self.args.ncore)
        # Map processing to _run function
        pool.imap(_run_skim, skim_files)
        # Close and join pool
        pool.close()
        pool.join()
