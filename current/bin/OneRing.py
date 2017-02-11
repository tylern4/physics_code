#!/usr/bin/env python
from __future__ import print_function
import argparse
import sys
import os
from ROOT import gBenchmark, gROOT
from pathos.multiprocessing import Pool, cpu_count, ProcessingPool
import multiprocessing as mp
try:
    import dill as pl
except:
    import cPickle as pl
else:
    import pickle as pl

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


def split_list(alist, wanted_parts=1):
    """found this on stack overflow"""
    length = len(alist)
    return [alist[i * length // wanted_parts: (i + 1) * length // wanted_parts]
            for i in range(wanted_parts)]


def skim(lines):
    for line in lines:
        infile = [line]
        outfile = line.replace('.root', '_skim.root').replace('/root', '/skim')
        command = "./Skim " + infile[0] + " " + outfile
        # print_yellow(command)
        os.system(command)


def main():
    parser = argparse.ArgumentParser(
        description="One Ring to Rule Them All and in the Darkness Bind Them:")
    parser.add_argument('input', type=str,
                        help="Input list for files")
    parser.add_argument('output', type=str, nargs='?',
                        help="Output folder for final root files", default='.')
    parser.add_argument('-n', dest='ncore', type=int, nargs='?',
                        help="Number of cores to use if not all the cores", default=0)
    # parser.add_argument('-r', dest='recompile', type=bool, nargs='?', action='store_true',
    #                    help="Recompile before starting run", default=False)

    args = parser.parse_args()
    if args.output[-1] != '/':
        args.output = args.output + '/'

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()

    os.system("cd ../cpp && make clean && make -j4 && cd ../bin")
    input_file = str(sys.argv[1])
    outfile_file = input_file[:-4]
    e1d_command = "./e1d " + input_file + " " + \
        outfile_file.replace('/inputFiles', '/outputFiles') + "_test.root"
    print_blue(e1d_command)
    os.system(e1d_command)

    print_red("e1d done")
    os.system("cd ../skim && make clean && make -j4 && cd ../bin")
    num_cores = (cpu_count(), args.ncore)[args.ncore > 0]
    pool = Pool(processes=num_cores)
    lines = [line.rstrip('\n') for line in open(input_file)]

    lines_split = split_list(lines, wanted_parts=num_cores)
    pool.map(skim, (lines_split))

    print_red("Skim done")

    input_file = str(sys.argv[1])
    outfile_file = input_file[:-4]
    e1d_command = "./e1d " + input_file[:-4] + "_skim.lis" + " " + \
        outfile_file.replace('/inputFiles', '/outputFiles') + "_cut.root 2"
    print_blue(e1d_command)
    os.system(e1d_command)
    print_red("Finished\a")

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print_red("\n\nExiting")
        sys.exit()
