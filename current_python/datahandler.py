import ROOT
import numpy as np 
from constants import *
from ROOT import TLorentzVector,TVector3
from ROOT import gBenchmark, gROOT
gROOT.SetBatch(True)
from physics import *
import plots
from pathos.multiprocessing import Pool, cpu_count, ProcessingPool
import multiprocessing as mp
try:
	import dill as pl
except:
	import cPickle as pl
else:
	import pickle as pl

import glob
import pandas as pd

from ROOT import gSystem
# load library with MyClass dictionary
gSystem.Load( 'kinematics_C' )
gSystem.Load( 'test_C' )
# get MyClass from ROOT
from ROOT import kinematics
from ROOT import test

class datahandeler(object):
	"""Datahandeler class"""
	def __init__(self, args):
		self.args = args
		self.Q2 = np.array([])
		self.W = np.array([])
		self.args.ncore = (cpu_count() , self.args.ncore)[self.args.ncore > 0]

	def split_list(self):
		wanted_parts = self.args.ncore
		alist = glob.glob(self.args.input + '*.root')
		length = len(alist)
		return [ alist[i*length // wanted_parts: (i+1)*length // wanted_parts]
				for i in range(wanted_parts) ]

	def _run(self, files):
		chain = ROOT.TChain('h10')
		t = test()
		for _f in files:
			chain.Add(_f)
		t.loop(chain)
		#del t


	def run_mp(self):
		files = self.split_list()
		pool = ProcessingPool()
		pool.map(self._run, (files))
		pool.close()
		pool.join()

		
	def plot(self):
		plot = plots.plotting()
		plot.WvsQ2(self.W,self.Q2, output=self.args.output)
