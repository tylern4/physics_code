from rootpy.tree import Tree
from rootpy.io import root_open
from ROOT import TChain, TTree, TBranch
import argparse
#from root_numpy import root2array, tree2array
import pandas as pd
import numpy as np
from tqdm import trange, trange
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

def getCurrentValue(event, branchName):
     return getattr(event,branchName)

#f = root_open(args.inList, "read")
parser = argparse.ArgumentParser(description='Test Root Funtionality')
parser.add_argument('inList', type=str, nargs='?')

args = parser.parse_args()

if args.inList is None:
	args.inList = "/Users/tylern/data/inputFiles/v3_skim.lis"

chain = TChain("h10") 

for line in open(args.inList,'r'): 
	chain.AddFile(line[:-1]) 

momentum = []
for event in trange(chain.GetEntries()):
	chain.GetEntry(event)
	p = getCurrentValue(chain, "p")
	for i in range(len(p)):
		momentum.append(p[i])

momentum = np.array(momentum)
n, bins, patches = plt.hist(momentum, 50, normed=1, edgecolor='None')
plt.show()