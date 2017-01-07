#!/usr/bin/env python
import ROOT
import sys
from main import *
from constants import *
from physics import *

file_name = '/Users/tylern/data/inputFiles/v3/skim/r22855_00_skim.root'
file_h10 = ROOT.TFile.Open(file_name, 'read')

#W = array
#Q2 = array

for _e in file_h10.h10:

	#append(W,_e.W)
	#append(Q2,_e.Q2)
	if _e.id[0] is ID['ELECTRON']:
		e_mu_p = fourvec(_e.p[0],_e.cx[0],_e.cy[0],_e.cz[0],mass['ELECTRON'])
		#Q2 = Q2_calc(e_mu,e_mu_p)
		print(_e.Q2 - Q2_calc(e_mu,e_mu_p))
	#for _i in xrange(len(_e.p)):
	#	px,py,pz = all_mom(_e.p[_i],_e.cx[_i],_e.cy[_i],_e.cz[_i])
	#	e_mu_p = fvec(px,py,pz,MASS_E)
	#	print((e_mu - e_mu_p).M())
		


