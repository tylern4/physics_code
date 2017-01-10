import ROOT
import numpy as np 
from constants import *
from ROOT import TLorentzVector,TVector3
from physics import *
import plots
from multiprocessing import Pool, cpu_count
import glob

class datahandeler(object):
	"""Datahandeler class"""
	def __init__(self, args,num_cores=0):
		self.args = args
		self.Q2 = np.array([])
		self.W = np.array([])
		self.num_cores = (cpu_count() , num_cores)[num_cores > 0]
		self.chain = ROOT.TChain('h10')
		self.list = self.split_list()

	def run(self):
		directory = self.args.input + '/*.root'
		self.chain.Add(directory)

		for _e in self.chain:
			if _e.id[0] is ID['ELECTRON']:
				e_mu_p = fourvec(_e.p[0],_e.cx[0],_e.cy[0],_e.cz[0],mass['ELECTRON'])
				self.Q2 = append(self.Q2,Q2_calc(e_mu,e_mu_p))
				self.W = append(self.W,W_calc(e_mu,e_mu_p))

	def plot(self):
		pl = plots.plotting()
		pl.WvsQ2(self.W,self.Q2, output=self.args.output)

	def split_list(self):
		wanted_parts = self.num_cores
		alist = glob.glob(self.args.input + '/*.root')
		length = len(alist)
		return [ alist[i*length // wanted_parts: (i+1)*length // wanted_parts] for i in range(wanted_parts) ]

