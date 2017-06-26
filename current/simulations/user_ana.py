#!/usr/bin/env python
from __future__ import print_function
from termcolor import cprint
print_green = lambda x: cprint(x, 'green', attrs=['bold'])
print_red = lambda x: cprint(x, 'red', attrs=['bold'])
print_blue = lambda x: cprint(x, 'blue', attrs=['bold'])
print_white = lambda x: cprint(x, 'white', attrs=['bold'])
print_yellow = lambda x: cprint(x, 'yellow', attrs=['bold'])

from multiprocessing import Pool, cpu_count
import multiprocessing as mp
import tqdm
import argparse
import sys
import glob
import contextlib
import os
import shutil
import tempfile

@contextlib.contextmanager
def cd(newdir, cleanup=lambda: True):
    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(newdir))
    try:
        yield
    finally:
        os.chdir(prevdir)
        cleanup()

@contextlib.contextmanager
def tempdir():
    dirpath = tempfile.mkdtemp(dir="/mnt/ssd/temp")
    def cleanup():
        shutil.rmtree(dirpath)
    with cd(dirpath, cleanup):
        yield dirpath

def split_list(args):
    """Split the list of input files into equal chunks to be processed by the pool."""
    wanted_parts = args.ncore
    alist = glob.glob(args.input + 'gsim*.evt')
    length = len(alist)
    return [alist[i * length // wanted_parts: (i + 1) * length // wanted_parts]
        for i in range(wanted_parts)]

def do_user_ana(base):
    with tempdir() as dirpath:
        shutil.copyfile(base, dirpath+"/uncooked.bos")
        shutil.copyfile("/home/tylern/physics_code/current/simulations/user_ana.tcl", dirpath+"/user_ana.tcl")
        shutil.copyfile("/home/tylern/physics_code/current/simulations/do_user_ana.sh", dirpath+"/do_user_ana.sh")
        os.system("docker run --link clasdb:clasdb -v`pwd`:/root/code --rm -it tylern4/clas6:latest do_user_ana.sh")
        shutil.copyfile(dirpath+"/cooked.root", base[:-4]+"_cooked.root")

def main():
    # Make argument parser
    parser = argparse.ArgumentParser(description="User analysis")
    parser.add_argument('input', type=str, help="Input directory for files to be user_ana'd")
    parser.add_argument('-n', dest='ncore', type=int, nargs='?', help="Number of cores to use if not all the cores", default=0)

    # Print help if there aren't enough arguments
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()

    args = parser.parse_args()
    # Make sure file paths have ending /
    if args.input[-1] != '/':
        args.input = args.input + '/'
    if args.ncore == 0 or cpu_count > cpu_count():
        args.ncore = cpu_count()

    #files = split_list(args)
    files = glob.glob(args.input + 'gsim*.evt')
    pool = Pool(processes=args.ncore)

    #output = pool.map(do_user_ana, files)
    for _ in tqdm.tqdm(pool.imap_unordered(do_user_ana, files), total=len(files)):
        pass

    # Close and join pool
    pool.close()
    pool.join()


if __name__ == "__main__":
    try:
	main()
    except KeyboardInterrupt:
        print_red("\n\nExiting")
        sys.exit()
