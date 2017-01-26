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
		kine = kinematics()
		t = test()
		for _f in files:
			chain.Add(_f)

		t.pass_chain(chain)
		#for _e in chain:
			#if _e.id[0] is ID['ELECTRON']:
			#	#e_mu_p = fourvec(_e.p[0],_e.cx[0],_e.cy[0],_e.cz[0],mass['ELECTRON'])
			#	kine.SetVal(_e.p[0],_e.cx[0],_e.cy[0],_e.cz[0])
			#	self.Q2 = append(self.Q2,kine.GetQ2())
			#	self.W = append(self.W,kine.GetW())
		#WQ2 = pd.DataFrame()
		#WQ2['W'] = self.W
		#WQ2['Q2'] = self.Q2

		#name = str(self.args.output+'WQ2_'+str(mp.current_process().pid)+'.pkl')
		#pl.dump(WQ2, open(name, "wb"),protocol=2)
		del kine
		del t


	def run_mp(self):
		files = self.split_list()
		pool = ProcessingPool()
		pool.map(self._run, (files))
		pool.close()
		pool.join()

		
	def plot(self):
		plot = plots.plotting()
		plot.WvsQ2(self.W,self.Q2, output=self.args.output)
