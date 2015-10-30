#!/usr/local/bin/ipython
from multiprocessing import Pool
import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from root_numpy import root2rec
from StringIO import StringIO
from ROOT import TLorentzVector
from ROOT import TVector3
from physics import Q2_calc, W_calc, branches, masses, mass, getM2
from missingmass import missing_mass_calc

def convert(lines):
	e_mu = TLorentzVector(0.0,0.0,np.sqrt(4.802**2.0 - 0.000511**2.0),4.802) #(0.0,0.0, sqrt(Square(E1D_E0)-Square(MASS_E)), E1D_E0)
	for line in lines:
		fileNames = [line]
		print fileNames
		
		df = pd.DataFrame(root2rec(fileNames, treename='h10', branches=branches)) #magic line to create a dataframe
	
		#These lines add a column to the data frame with 'name' = value, should be able to find W, Q2 for an event and add it like this
		df['pz'] = df.p*df.cz 
		df['py'] = df.p*df.cy
		df['px'] = df.p*df.cx
		df['energy'] = (df.p*df.p + df.id.apply(getM2)).apply(lambda x: np.sqrt(x)) #This is some master python
		
		#more complicated calcs and if statement ones can be added like this
		Q2 = []
		W = []
		MM = []
		np.array(Q2)
		np.array(W)
		np.array(MM)
		gamma_mu = TLorentzVector()
		pi_mu = TLorentzVector()
		for i in xrange(len(df)):
			#df['stat'][0] > 0 and df['q'][0] == -1 and df['sc'][0] > 0 and df['dc'][0] > 0 and df['ec'][0] > 0 and df['dc_stat'][df['dc'][0]-1] > 0:
			if df['id'][i][0] == 11 and df['gpart'][i] > 1:
				e_mu_prime_3 = TVector3(df.px[i][0],df.py[i][0],df.pz[i][0])	
				e_mu_prime = TLorentzVector()
				e_mu_prime.SetVectM(e_mu_prime_3, mass['ELECTRON'])
				Q2.append(Q2_calc(e_mu, e_mu_prime))
				W.append(W_calc(e_mu,e_mu_prime))
	
				if df['id'][i][1] == 211:
					gamma_mu = (e_mu - e_mu_prime)
					pi_mu_3 = TVector3(df.px[i][1],df.py[i][1],df.pz[i][1])
					pi_mu.SetVectM(pi_mu_3, mass['PIP'])
					mm = missing_mass_calc(gamma_mu, pi_mu)
					MM.append(mm)
				else:
					MM.append(np.NaN)
			else:
				Q2.append(np.NaN)
				W.append(np.NaN)
				MM.append(np.NaN)
	
		df['Q2'] = pd.Series(Q2, index=df.index)
		df['W'] = pd.Series(W, index=df.index)
		df['MM'] = pd.Series(MM, index=df.index)
	
		df = df[~(df.Q2 != df.Q2)]
	
	
	
		fileNames = line.replace('root', 'h5')
		os.remove(fileNames) if os.path.isfile(fileNames) else None
		store = pd.HDFStore(fileNames)
		store.put('df',df)
		store.close()


lines = [line.rstrip('\n') for line in open(str(sys.argv[1]))]
num_cores = 8
pool = Pool(processes=num_cores)
iteration = len(lines)/num_cores

lines1 = lines[:iteration]
lines2 = lines[1*iteration+1:2*iteration]
lines3 = lines[2*iteration+1:3*iteration]
lines4 = lines[3*iteration+1:4*iteration]
lines5 = lines[4*iteration+1:5*iteration]
lines6 = lines[5*iteration+1:6*iteration]
lines7 = lines[6*iteration+1:7*iteration]
lines8 = lines[7*iteration+1:]
length = len(lines1) + len(lines2) + len(lines3) + len(lines4) + len(lines5) + len(lines6) + len(lines7) + len(lines8)
if length == len(lines):
	pool.map(convert, (lines1,lines2,lines3,lines4,lines5,lines6,lines7,lines8))
else:
	print "error: length ", length, "!=", len(lines) 

print "done"