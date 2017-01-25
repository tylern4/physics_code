import ROOT
import numpy as np 
from constants import *
from ROOT import TLorentzVector,TVector3
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
#from multiprocessing import Pool, cpu_count
import glob
from ROOT import gBenchmark
import pandas as pd

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

	def run(self):
		chain = ROOT.TChain('h10')
		fnum = chain.Add(self.args.input + '*.root')
		print("Processing "+ str(chain.GetEntries())+ " events in "+ str(fnum) + " files")
		for _e in chain:
			if _e.id[0] is ID['ELECTRON']:
				e_mu_p = fourvec(_e.p[0],_e.cx[0],_e.cy[0],_e.cz[0],mass['ELECTRON'])
				self.Q2 = append(self.Q2,Q2_calc(e_mu,e_mu_p))
				self.W = append(self.W,W_calc(e_mu,e_mu_p))

	def _run(self, files):
		chain = ROOT.TChain('h10')
		for _f in files:
			chain.Add(_f)
		for _e in chain:
			if _e.id[0] is ID['ELECTRON']:
				e_mu_p = fourvec(_e.p[0],_e.cx[0],_e.cy[0],_e.cz[0],mass['ELECTRON'])
				self.Q2 = append(self.Q2,Q2_calc(e_mu,e_mu_p))
				self.W = append(self.W,W_calc(e_mu,e_mu_p))
		WQ2 = pd.DataFrame()
		WQ2['W'] = self.W
		WQ2['Q2'] = self.Q2

		name = str(self.args.output+'WQ2_'+str(mp.current_process().pid)+'.pkl')
		pl.dump(WQ2, open(name, "wb"),protocol=2)


	def run_mp(self):
		files = self.split_list()
		pool = ProcessingPool()
		pool.map(self._run, (files))
		pool.close()
		pool.join()

		
	def plot(self):
		plot = plots.plotting()
		plot.WvsQ2(self.W,self.Q2, output=self.args.output)
