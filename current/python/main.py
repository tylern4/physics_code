#!/usr/bin/env python
from __future__ import print_function
import argparse
import sys
from ROOT import gBenchmark, gROOT
from datahandler import datahandeler
from constants import *

gROOT.SetBatch(True)


def main():
    # Make argument parser
    parser = argparse.ArgumentParser(description="Root datahandeler program")
    parser.add_argument('input', type=str,
                        help="Input directory for *.root files")
    parser.add_argument('output', type=str, nargs='?',
                        help="Output for pdf files", default='.')
    parser.add_argument('-n', dest='ncore', type=int, nargs='?',
                        help="Number of cores to use if not all the cores", default=0)
    parser.add_argument('-c', dest='clean',
                        help="Clean output directory of all pkl before start", default=False, action='store_true')

    # Print help if there aren't enough arguments
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()

    args = parser.parse_args()

    # Make sure file paths have ending /
    if args.input[-1] != '/':
        args.input = args.input + '/'
    if args.output[-1] != '/':
        args.output = args.output + '/'

    # Create datahandeler object based on the given arguments
    dh = datahandeler(args)

    gBenchmark.Start('Map')
    dh.run_map()
    gBenchmark.Show('Map')

    # TODO: Get plotting function to work
    gBenchmark.Start('Reduce')
    dh.run_reduce()
    gBenchmark.Show('Reduce')


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print_red("\n\nExiting")
        sys.exit()
