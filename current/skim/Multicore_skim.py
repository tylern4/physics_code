#!/usr/bin/env python
from skim import skim, split_list
from multiprocessing import Pool
import multiprocessing
import sys
import os
from timeit import default_timer as timer

try:
    start = timer()
    os.system("make clean && make")

    lines = [line.rstrip('\n') for line in open(str(sys.argv[1]))]

    num_cores = multiprocessing.cpu_count()
    pool = Pool(processes=num_cores)

    lines_split = split_list(lines, wanted_parts=num_cores)

    pool.map(skim, (lines_split))
    end = timer()
    print(end - start)
except KeyboardInterrupt:
    print("Why do you want to kill me??")
    pool.terminate()
    pool.join()
    sys.exit(0)


print("done")
