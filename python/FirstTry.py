#!/usr/bin/env python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from root_numpy import root2rec, root2array
from StringIO import StringIO
from ROOT import TLorentzVector

masses = {11:0.000511, 211:0.13957, -211:0.13957, 2212:0.93827, 2112:0.939565, 0:0, 22:0, 321:0.493667, -321:0.493667, 45:0, 47:0, 49:0}
getM2 = lambda x: [masses[pid]**2 for pid in x]

branches = ['npart',
			'gpart', 
			'id',
			'evntid',
			'evtype',
			'evntclas',
			'evthel',
			'evntclas2',
			'q_l',
			't_l',
			'tr_time',
			'rf_time1',
			'rf_time2',
			
			'stat',   
			'dc',   
			'cc',   
			'sc',   
			'ec',   
			'lec',   
			'p',   
			'q',   
			'b',   
			'cx',   
			'cy',   
			'cz',   
			'vx',   
			'vy',   
			'vz',   
			
			'dc_part',		
			'dc_sect',   
			'dc_trk',   
			'dc_stat',   
			'dc_xsc',   
			'dc_ysc',   
			'dc_zsc',   
			'dc_cxsc',   
			'dc_cysc',   
			'dc_czsc',   
			'dc_xec',   
			'dc_yec',   
			'dc_zec',   
			'dc_thcc',   
			'dc_c2',   
			
			'ec_part',
			'ec_stat',   
			'ec_sect',   
			'ec_whol',   
			'etot',   
			'ec_ei',   
			'ec_eo',   
			'ec_t',   
			'ec_r',   
			'ech_x',   
			'ech_y',   
			'ech_z',   
			'ec_c2',   
			
			'sc_part',
			'sc_sect',   
			'sc_hit',   
			'sc_pd',   
			'sc_stat',   
			'edep',   
			'sc_t',   
			'sc_r',   
			'sc_c2',   
			
			'cc_part',
			'cc_sect',   
			'cc_hit',   
			'cc_segm',   
			'nphe',   
			'cc_t',   
			'cc_r',   
			'cc_c2'
			] #This sets the branch names, should be able to add all of them to load everything or just use the ones I need for now

#n, bins, patches = plt.hist(energy, 500, normed=0, range=(0.001, 3.14), histtype=u'stepfilled', facecolor='g' , alpha=0.45)
#plt.show()
lines = [line.rstrip('\n') for line in open('v3_all.lis')]

for line in lines:
	fileNames = [line]
	print fileNames
	df = pd.DataFrame(root2rec(fileNames, treename='h10', branches=branches)) #magic line to create a dataframe
	#These lines add a column to the data frame with 'name' = value, should be able to find W, Q2 for an event and add it like this
	df['pz'] = df.p*df.cz 
	df['py'] = df.p*df.cy
	df['px'] = df.p*df.cx
	df['energy'] = (df.p*df.p + df.id.apply(getM2)).apply(lambda x: np.sqrt(x)) #This is some master python





	fileNames = line.replace('root', 'h5')
	store = pd.HDFStore(fileNames)
	store.put('df',df)

#energy = np.hstack(np.array(df['energy'])) ###Flattens array or arrays to a single array