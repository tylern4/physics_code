#!/usr/local/bin/ipython
import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from root_numpy import root2rec
from StringIO import StringIO
from physics import branches

lines = [line.rstrip('\n') for line in open(str(sys.argv[1]))]

for line in lines:
	fileNames = [line]
	print fileNames
	
	df = pd.DataFrame(root2rec(fileNames, treename='h10', branches=branches))
	print df