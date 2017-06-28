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

def make_list(num):
    from datetime import datetime
    time = datetime.now().strftime('%m_%d_%Y-%H%M_')
    l = []
    for i in range(0,num):
        l.append("sim_"+time+str(i))
    return l

def do_sim(base):
    with tempdir() as dirpath:
        shutil.copyfile("/home/tylern/physics_code/current/simulations/aao_rad.inp", dirpath+"/aao_rad.inp")
        shutil.copyfile("/home/tylern/physics_code/current/simulations/gsim.inp", dirpath+"/gsim.inp")
        shutil.copyfile("/home/tylern/physics_code/current/simulations/user_ana.tcl", dirpath+"/user_ana.tcl")

        shutil.copyfile("/home/tylern/physics_code/current/simulations/do_sim.sh", dirpath+"/do_sim.sh")

        os.system("docker run --link clasdb:clasdb -v`pwd`:/root/code --rm -it tylern4/clas6:latest do_sim.sh")
        shutil.copyfile(dirpath+"/cooked.root", base+".root")

def main():
    # Make argument parser
    parser = argparse.ArgumentParser(description="Full sim analysis")
    parser.add_argument('output', type=str, help="Output directory for final root files")
    parser.add_argument('-n', dest='ncore', type=int, nargs='?', help="Number of cores to use if not all the cores", default=0)
    parser.add_argument('-c', dest='count', type=int, nargs='?', help="Number of passes to do", default=0)

    # Print help if there aren't enough arguments
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()

    args = parser.parse_args()
    # Make sure file paths have ending /
    if args.output[-1] != '/':
        args.output = args.output + '/'
    if args.ncore == 0 or cpu_count > cpu_count():
        args.ncore = cpu_count()

    files = make_list(args.count)
    pool = Pool(processes=args.ncore)
    print(files)
    #output = pool.map(do_user_ana, files)
    for _ in tqdm.tqdm(pool.imap_unordered(do_sim, files), total=len(files)):
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
