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
gSystem.Load( 'H10_C' )
# get MyClass from ROOT
from ROOT import H10
from ROOT import std
import array

class datahandeler(object):
	"""Datahandeler class"""
	def __init__(self, args):
		self.args = args
		self.Q2 = np.array([])
		self.W = np.array([])
		self.W_Q2 = pd.DataFrame()
		self.p_beta = pd.DataFrame()
		self.args.ncore = (cpu_count() , self.args.ncore)[self.args.ncore > 0]
		print("Starting datahandeler with %d cores" % self.args.ncore)

	def split_list(self):
		wanted_parts = self.args.ncore
		alist = glob.glob(self.args.input + '*.root')
		length = len(alist)
		return [ alist[i*length // wanted_parts: (i+1)*length // wanted_parts]
				for i in range(wanted_parts) ]

	def _run(self, files):
		chain = ROOT.TChain('h10')
		chain.UseCache(100,1024)
		h10 = H10()
		for _f in files:
			chain.Add(_f)
		h10.loop(chain)
		print("Done with loop "+str(mp.current_process().pid))
		#_p_beta['p'] = [_p for _p in h10.p_vec]
		#_p_beta['b'] = [_b for _b in h10.b_vec]
		#_p_beta['q'] = [_q for _q in h10.q_vec]
		#_p_beta['id'] = [_id for _id in h10.id_vec]
		#pl.dump(_p_beta, open(self.args.output + 'p_beta_'+str(mp.current_process().pid)+'.pkl','wb'))
		#del _p_beta 
		#print("Done with _p_beta "+str(mp.current_process().pid))
		#_W_Q2['W'] = [_W for _W in h10.W_vec]
		#_W_Q2['Q2'] = [_Q2 for _Q2 in h10.Q2_vec]
		#pl.dump(_W_Q2, open(self.args.output + 'W_Q2_'+str(mp.current_process().pid)+'.pkl','wb'))
		#del _W_Q2
		#print("Done with _W_Q2 "+str(mp.current_process().pid))
		#del h10
		#print("Done with "+str(mp.current_process().pid))
		return h10



	def run_mp(self):
		files = self.split_list()
		pool = ProcessingPool()
		out = pool.map(self._run, (files))
		pool.close()
		pool.join()
		i = 0
		for _h10 in out:
			i = i+1
			print("Doing loop "+str(i))
			_p_beta = pd.DataFrame()
			_p_beta['p'] = [_p for _p in _h10.p_vec]
			_p_beta['b'] = [_b for _b in _h10.b_vec]
			_p_beta['q'] = [_q for _q in _h10.q_vec]
			_p_beta['id'] = [_id for _id in _h10.id_vec]
			pl.dump(_p_beta, open(self.args.output + 'p_beta_'+str(i)+'.pkl','wb'))
			del _p_beta
			_W_Q2 = pd.DataFrame()
			_W_Q2['W'] = [_W for _W in _h10.W_vec]
			_W_Q2['Q2'] = [_Q2 for _Q2 in _h10.Q2_vec]
			pl.dump(_W_Q2, open(self.args.output + 'W_Q2_'+str(i)+'.pkl','wb'))
			del _W_Q2

		#self.W_Q2 = pd.concat(out, ignore_index=False)

	def plot(self):
		plot = plots.plotting()
		plot.WvsQ2(self.W_Q2, output=self.args.output)
