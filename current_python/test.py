from rootpy.tree import Tree
from rootpy.io import root_open
from ROOT import TChain, TTree, TBranch
import argparse

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


for event in range(chain.GetEntries()):
	chain.GetEntry(event)
	p = getCurrentValue(chain, "p")
	for i in range(len(p)):
		print(p[i])
