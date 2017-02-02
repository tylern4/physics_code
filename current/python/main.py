#!/usr/bin/env python
from __future__ import print_function
from datahandler import datahandeler
import argparse
import sys
from ROOT import gBenchmark, gROOT

try:
    from termcolor import cprint
    print_green = lambda x: cprint(x, 'green', attrs=['bold'])
    print_red = lambda x: cprint(x, 'red', attrs=['bold'])
    print_blue = lambda x: cprint(x, 'blue', attrs=['bold'])
    print_white = lambda x: cprint(x, 'white', attrs=['bold'])
    print_yellow = lambda x: cprint(x, 'yellow', attrs=['bold'])
except ImportError:
    print_green = lambda x: print(x)
    print_red = lambda x: print(x)
    print_blue = lambda x: print(x)
    print_white = lambda x: print(x)
    print_yellow = lambda x: print(x)

gROOT.SetBatch(True)


def main():
    parser = argparse.ArgumentParser(description="Root datahandeler program")
    parser.add_argument('input', type=str,
                        help="Input directory for *.root files")
    parser.add_argument('output', type=str, nargs='?',
                        help="Output for pdf files", default='.')
    parser.add_argument('-n', dest='ncore', type=int, nargs='?',
                        help="Number of cores to use if not all the cores", default=0)

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()

    args = parser.parse_args()
    if args.input[-1] != '/':
        args.input = args.input + '/'
    if args.output[-1] != '/':
        args.output = args.output + '/'

    dh = datahandeler(args)

    gBenchmark.Start('Run')
    dh.run_mp()
    gBenchmark.Show('Run')

    # gBenchmark.Start('Plot')
    # dh.plot()
    # gBenchmark.Show('Plot')


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print_red("\n\nExiting")
        sys.exit()
