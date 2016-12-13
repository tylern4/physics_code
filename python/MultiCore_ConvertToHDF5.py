#!/usr/bin/env python
from convert import convert, split_list
from multiprocessing import Pool
import multiprocessing
import sys

lines = [line.rstrip('\n') for line in open(str(sys.argv[1]))]

num_cores = multiprocessing.cpu_count()
pool = Pool(processes=num_cores)

lines_split = split_list(lines, wanted_parts=num_cores)

pool.map(convert, (lines_split))

print("done")
