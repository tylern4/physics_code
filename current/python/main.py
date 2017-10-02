#!/usr/bin/env python
from __future__ import print_function, unicode_literals, absolute_import
import argparse
import sys
from ROOT import gBenchmark, gROOT
from datahandler import datahandeler
from constants import *

gROOT.SetBatch(True)
from halo import Halo


def main():
    spinner = Halo({'text': 'Running', 'spinner': 'dots'})
    # Make argument parser
    parser = argparse.ArgumentParser(description="Root datahandeler program")
    parser.add_argument('input', type=str,
                        help="Input directory for *.root files")
    parser.add_argument('output', type=str, nargs='?',
                        help="Output root file", default='./test.root')
    parser.add_argument('-n', dest='ncore', type=int, nargs='?',
                        help="Number of cores to use if not all the cores", default=0)
    parser.add_argument('-c', dest='clean',
                        help="Clean output directory of all pkl before start", default=False, action='store_true')
    parser.add_argument('-s', dest='skim',
                        help="Also skim files and run analysis", default=False, action='store_true')

    # Print help if there aren't enough arguments
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()

    args = parser.parse_args()
    spinner.start()
    # Make sure file paths have ending /
    if args.input[-1] != '/':
        args.input = args.input + '/'
    if args.output[-5:] != '.root':
        args.output = args.output + '.root'

    spinner.succeed('Starting DataHandeler')
    spinner.start()
    # Create datahandeler object based on the given arguments
    dh = datahandeler(args)

    spinner.succeed('Starting Map')
    spinner.start()
    gBenchmark.Start('Map')
    dh.run_map()
    gBenchmark.Show('Map')

    spinner.succeed('Starting Reduce')
    spinner.start()
    gBenchmark.Start('Reduce')
    dh.run_reduce()
    gBenchmark.Show('Reduce')
    spinner.stop()

    # TODO get skim to work properly
    if args.skim:
        gBenchmark.Start('Skim')
        dh.run_skim()
        gBenchmark.Show('Skim')

        #args.input = args.input.replace("/root", "/skim")
        #args.output = args.output.replace(".root", "_skim.root")
        #skim = datahandeler(args)

        # gBenchmark.Start('Map_Skim')
        # skim.run_map()
        # gBenchmark.Show('Map_Skim')

        # gBenchmark.Start('Reduce_Skim')
        # skim.run_reduce()
        # gBenchmark.Show('Reduce_Skim')


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print_red("\n\nExiting")
        sys.exit()
