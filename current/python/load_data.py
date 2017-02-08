#!/usr/bin/env python
import numpy as np
import pandas as pd
try:
    import dill as pl
except:
    import cPickle as pl
else:
    import pickle as pl
import argparse
import glob
import sys
import plots
from ROOT import gBenchmark, gROOT
gROOT.SetBatch(True)


def load_values():
    WQ2 = []
    all_WQ2_pickles = glob.glob(args.input + 'WQ2_*.pkl')
    for _file in all_WQ2_pickles:
        with open(_file, 'rb') as in_strm:
            WQ2.append(pl.load(in_strm))

    return pd.concat(WQ2, ignore_index=False)


def load_graphs():
    with open(args.output + 'WvsQ2.pdf.pkl', 'rb') as in_strm:
        fig = pl.load(in_strm)
    plt.show()


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Root datahandeler program")
    parser.add_argument('input', type=str,
                        help="Input directory for *.pkl files")
    parser.add_argument('output', type=str, nargs='?',
                        help="Output for pdf files", default='.')
    parser.add_argument('-s', help="Save outputs to pdf",
                        default=False, action='store_true')
    parser.add_argument('-p', help="Save outputs to pickle",
                        default=False, action='store_true')
    parser.add_argument('-o', help="Open outputs as plots",
                        default=False, action='store_true')

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()

    args = parser.parse_args()
    if args.input[-1] != '/':
        args.input = args.input + '/'
    if args.output[-1] != '/':
        args.output = args.output + '/'

    if args.s or args.p:
        WQ2 = load_values()
        plot = plots.plotting(save=args.s, pickle=args.p)
        plot.WvsQ2(WQ2['W'], WQ2['Q2'], output=args.output)

    if args.o:
        fig = load_graphs()
        plt.show()
