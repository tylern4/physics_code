import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from root_numpy import root2rec
from ROOT import TLorentzVector
from ROOT import TVector3
from physics import Q2_calc, W_calc, branches, masses, mass, getM2
from missingmass import missing_mass_calc
import xray

#found this on stack overflow
def split_list(alist, wanted_parts=1):
    length = len(alist)
    return [ alist[i*length // wanted_parts: (i+1)*length // wanted_parts] 
             for i in range(wanted_parts) ]

def convert(lines):
	for line in lines:

		fileName = [line]
		print(fileName)
		#magic line to create a dataframe from ROOT file
		#df = pd.DataFrame(root2array(fileName, treename='h10', branches=branches).view(np.recarray))
		df = pd.DataFrame(root2rec(fileName, treename='h10', branches=branches)) 
		df = df[df.W >= 0]

		fileName = line.replace('root', 'h5').replace('/skim','/h5')
		os.remove(fileName) if os.path.isfile(fileName) else None
		store = pd.HDFStore(fileName,format='table')
		store.put('df', df)
		store.close()