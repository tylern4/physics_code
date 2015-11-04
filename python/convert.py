import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from root_numpy import root2rec
from ROOT import TLorentzVector
from ROOT import TVector3
from physics import Q2_calc, W_calc, branches, masses, mass, getM2
from missingmass import missing_mass_calc

#found this on stack overflow
def split_list(alist, wanted_parts=1):
    length = len(alist)
    return [ alist[i*length // wanted_parts: (i+1)*length // wanted_parts] 
             for i in range(wanted_parts) ]

def convert(lines):
	e_mu = TLorentzVector(0.0,0.0,np.sqrt(4.802**2.0 - 0.000511**2.0),4.802) #(0.0,0.0, sqrt(Square(E1D_E0)-Square(MASS_E)), E1D_E0)
	for line in lines:
		fileNames = [line]

		print fileNames
		df = pd.DataFrame(root2rec(fileNames, treename='h10', branches=branches)) #magic line to create a dataframe
		
		#These lines add a column to the data frame with 'name' = value, should be able to find W, Q2 for an event and add it like this
		df['px'] = df.p*df.cx
		df['py'] = df.p*df.cy
		df['pz'] = df.p*df.cz 
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