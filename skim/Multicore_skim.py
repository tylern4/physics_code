#!/usr/local/bin/ipython
from skim import skim, split_list
# from StringIO import StringIO
from multiprocessing import Pool
import multiprocessing
import sys
import os

os.system("make")

lines = [line.rstrip('\n') for line in open(str(sys.argv[1]))]

num_cores = multiprocessing.cpu_count()
pool = Pool(processes=num_cores)

lines_split = split_list(lines, wanted_parts=num_cores)

pool.map(skim, (lines_split))

print "done"
