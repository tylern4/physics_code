#!/usr/bin/env python
from __future__ import print_function
import argparse
import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
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


def main():
    parser = argparse.ArgumentParser(
        description="Plot data from golden run output csv:")
    parser.add_argument('input', type=str,
                        help="Input CSV")
    parser.add_argument('output', type=str, nargs='?',
                        help="Output folder for final root files", default='.')

    args = parser.parse_args()
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()
    if args.input[-3:] != 'csv':
        print_white(args.input[-3])
        print_red("Error: Input must be CSV file")
        parser.print_help()
        parser.exit()
    if args.output[-1] != '/':
        args.output = args.output + '/'

    golden_df = pd.read_csv(args.input)
    golden_df['ratio'] = golden_df['num_of_events'] / golden_df['total_q']

    fig = plt.figure(num=None, figsize=(16, 9), dpi=200,
                     facecolor='w', edgecolor='k')
    plt.scatter(golden_df['run_num'], golden_df['ratio'], alpha=0.1)
    plt.ylim([-10000, 150000])
    plt.savefig('golden_run.pdf')


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print_red("\n\nExiting")
        sys.exit()