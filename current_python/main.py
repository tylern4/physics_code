#!/usr/bin/env python
from datahandler import datahandeler
import argparse
import sys
from ROOT import gBenchmark

def main():
	parser = argparse.ArgumentParser(description="Root datahandeler program")
	parser.add_argument('input', type=str, help="Input directory for *.root files")
	parser.add_argument('output', type=str, nargs='?', help="Output for pdf files", default='.')

	if len(sys.argv[1:])==0:
		parser.print_help()
		parser.exit()
	
	gBenchmark.Start('Args')
	args = parser.parse_args()

	#dh = datahandeler(args, num_cores=1)
	dh = datahandeler(args)
	gBenchmark.Show('Args')

	gBenchmark.Start('Run')
	dh.run()
	gBenchmark.Show('Run')

	gBenchmark.Start('Plot')
	dh.plot()
	gBenchmark.Show('Plot')



if __name__ == "__main__": 
	main()