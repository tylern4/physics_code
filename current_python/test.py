import ROOT
from ROOT import gROOT
from array import array
from physics import *
import numpy as np
from ROOT import gBenchmark

gROOT.SetBatch(True)
chain = ROOT.TChain('h10')

def defineArray(type,size):
	obj = array(type)
	for f in range(0,size):
		obj.append(0)
	return obj

chain.Add('/Users/tylern/Desktop/dumb/*.root')

#get the entries
entries = chain.GetEntries()

#declaration of the variables
p = array( 'd', [0] )
chain.SetBranchAddress("p",p)
cx = array( 'd', [0] )
chain.SetBranchAddress("cx",cx)
cy = array( 'd', [0] )
chain.SetBranchAddress("cy",cy)
cz = array( 'd', [0] )
chain.SetBranchAddress("cz",cz)
#ids = array( 'I', [0] )
#chain.SetBranchAddress("id",ids)
chain.SetBranchStatus('*',1)

Q2 = np.array([])
W = np.array([])

#getting the number of events
nEvt=chain.GetEntriesFast()
print(nEvt)
#looping over events
gBenchmark.Start('Run')
for jentry in xrange(0,nEvt):
	chain.GetEntry(jentry)
	#if id[0] is ID['ELECTRON']:
		#print("Hi")
	e_mu_p = fourvec(p[0],cx[0],cy[0],cz[0],mass['ELECTRON'])
	Q2 = append(Q2,Q2_calc(e_mu,e_mu_p))
	W = append(W,W_calc(e_mu,e_mu_p))
	#print(jentry)
gBenchmark.Show('Run')
#print(W,Q2)
#del chain
