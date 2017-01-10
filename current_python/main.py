#!/usr/bin/env python
from datahandler import datahandeler
import argparse
import sys

def main():
	parser = argparse.ArgumentParser(description="Root datahandeler program")
	parser.add_argument('input', type=str, help="Input directory for *.root files")
	parser.add_argument('output', type=str, nargs='?', help="Output for pdf files", default='.')

	if len(sys.argv[1:])==0:
		parser.print_help()
		parser.exit()

	args = parser.parse_args()

	dh = datahandeler(args)
	dh.run()
	dh.plot()


if __name__ == "__main__": 
	main()