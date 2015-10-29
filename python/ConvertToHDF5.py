#!/usr/local/bin/ipython
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from root_numpy import root2rec
from StringIO import StringIO
from ROOT import TLorentzVector
from ROOT import TVector3
from physics import Q2_calc, W_calc, branches, masses, getM2

lines = [line.rstrip('\n') for line in open(str(sys.argv[1]))]
e_mu = TLorentzVector(0.0,0.0,4.802,4.802) #(0.0,0.0, sqrt(Square(E1D_E0)-Square(MASS_E)), E1D_E0)
for line in lines:
	fileNames = [line]
	print fileNames
	
	df = pd.DataFrame(root2rec(fileNames, treename='h10', branches=branches)) #magic line to create a dataframe
	#These lines add a column to the data frame with 'name' = value, should be able to find W, Q2 for an event and add it like this
	df['pz'] = df.p*df.cz 
	df['py'] = df.p*df.cy
	df['px'] = df.p*df.cx
	df['energy'] = (df.p*df.p + df.id.apply(getM2)).apply(lambda x: np.sqrt(x)) #This is some master python
	

	Q2 = []
	W = []
	np.array(Q2)
	np.array(W)
	for i in xrange(len(df)):
		#df['stat'][0] > 0 and df['q'][0] == -1 and df['sc'][0] > 0 and df['dc'][0] > 0 and df['ec'][0] > 0 and df['dc_stat'][df['dc'][0]-1] > 0:
		if df['id'][i][0] == 11 and df['gpart'][i] > 0:
			e_mu_prime_3 = TVector3(df.px[i][0],df.py[i][0],df.pz[i][0])	
			e_mu_prime = TLorentzVector()
			e_mu_prime.SetVectM(e_mu_prime_3, 0.000511)
			Q2.append(Q2_calc(e_mu, e_mu_prime))
			W.append(W_calc(e_mu,e_mu_prime))
		else:
			Q2.append(-1)
			W.append(-1)
	df['Q2'] = pd.Series(Q2, index=df.index)
	df['W'] = pd.Series(W, index=df.index)


	fileNames = line.replace('root', 'h5')
	store = pd.HDFStore(fileNames)
	store.put('df',df)
	store.close()

#energy = np.hstack(np.array(df['energy'])) ###Flattens array or arrays to a single array
#Q2 = np.hstack(np.array(df['Q2']))
#W = np.hstack(np.array(df['W']))
#Q2[Q2 >= 0]
#W[W >= 0]
#plt.hist2d(W, Q2,bins=500,range=[[0,3.14],[0,4]])
#plt.show()
