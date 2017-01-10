import ROOT
import numpy as np 
from constants import *
from ROOT import TLorentzVector,TVector3
from physics import *
import plots

class datahandeler(object):
	"""Datahandeler class"""
	def __init__(self, args):
		self.args = args
		self.Q2 = np.array([])
		self.W = np.array([])

	def run(self):
		directory = self.args.input + '/*.root'
		chain_h10 = ROOT.TChain('h10')
		chain_h10.Add(directory)

		for _e in chain_h10:
			if _e.id[0] is ID['ELECTRON']:
				e_mu_p = fourvec(_e.p[0],_e.cx[0],_e.cy[0],_e.cz[0],mass['ELECTRON'])
				self.Q2 = append(self.Q2,Q2_calc(e_mu,e_mu_p))
				self.W = append(self.W,W_calc(e_mu,e_mu_p))

		
	def plot(self):
		pl = plots.plotting()
		pl.WvsQ2(self.W,self.Q2, output=self.args.output)