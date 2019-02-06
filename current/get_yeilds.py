#!/usr/bin/env python
import numpy as np

import glob
import sys
import time
from multiprocessing import Pool, cpu_count
from tqdm import tqdm
import argparse

import yeilds


def _yeilds(input_root):
    output_csv = input_root + ".csv"
    input_root = input_root + "*.root"
    yeilds.yeild(input_root.encode("utf-8"), output_csv.encode("utf-8")).run()


def main(args):
    if "*" in args.input:
        files = glob.glob(args.input)
    elif "/" in args.input:
        files = glob.glob(args.input + "/*.root")

    base = files[0][: files[0].rfind("/") + 1]
    files = [f[f.rfind("/") + 1 :] for f in files]
    files = [f[: f.rfind(".A")] for f in files]
    files = [f[: f.rfind("_")] for f in files]
    files = set(files)
    files = [base + f for f in files]

    time_start = time.time()
    num_cores = cpu_count()
    pool = Pool(processes=num_cores)
    _proc = pool.map(_yeilds, (files))
    pool.close()
    time_end = time.time()
    print("Time: {} Sec".format(time_end - time_start))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Python ROOT program")
    parser.add_argument("input", type=str, help="Input directory for *.root files")

    # Print help if there aren't enough arguments
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()

    args = parser.parse_args()
    main(args)
