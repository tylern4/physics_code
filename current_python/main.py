#!/usr/bin/env python
from datahandler import datahandeler
import argparse
import sys
from ROOT import gBenchmark, gROOT
gROOT.SetBatch(True)

def main():
	parser = argparse.ArgumentParser(description="Root datahandeler program")
	parser.add_argument('input', type=str, help="Input directory for *.root files")
	parser.add_argument('output', type=str, nargs='?', help="Output for pdf files", default='.')
	parser.add_argument('-n', dest='ncore', type=int, nargs='?', help="Number of cores to use if not all the cores", default=0)

	if len(sys.argv[1:])==0:
		parser.print_help()
		parser.exit()
	
	args = parser.parse_args()
	if args.input[-1] != '/':
		args.input = args.input+'/'
	if args.output[-1] != '/':
		args.output = args.output+'/'

	dh = datahandeler(args)

	gBenchmark.Start('Run')
	dh.run_mp()
	gBenchmark.Show('Run')

	gBenchmark.Start('Plot')
	dh.plot()
	gBenchmark.Show('Plot')



if __name__ == "__main__":
	try:
		main()
	except KeyboardInterrupt:
		print("Exiting")
		sys.exit()
