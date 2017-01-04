#!/usr/bin/env python
from __future__ import print_function
import argparse
from multiprocessing import Pool
import multiprocessing
import os
import sys

#found this on stack overflow
def split_list(alist, wanted_parts=1):
    length = len(alist)
    return [ alist[i*length // wanted_parts: (i+1)*length // wanted_parts] 
             for i in range(wanted_parts) ]

def h2root(lines):
	for line in lines:
		infile = [line]
		outfile = line.replace('hbook', 'root')
		command = "h2root "+infile[0]+" "+outfile
		os.system(command)
try:
	parser = argparse.ArgumentParser(description='Multicore h2root')
	parser.add_argument('input', help="Input Directory", nargs='?',type=str)
	
	args = parser.parse_args()
	
	if args.input is None:
		args.input = '.'
	
	files = [f for f in os.listdir(args.input) if os.path.isfile(args.input+f)]
	lines = [args.input+f for f in files if "hbook" in f]
	
	num_cores = multiprocessing.cpu_count()
	pool = Pool(processes=num_cores)
	lines_split = split_list(lines, wanted_parts=num_cores)
	pool.map(h2root, (lines_split))

except:
	e = sys.exc_info()[0]
	print(e)
