#!/usr/local/bin/ipython
from multiprocessing import Pool
import multiprocessing
import sys
import os
from termcolor import colored,cprint
print_green = lambda x: cprint(x,'green',attrs=['bold'])
print_red = lambda x: cprint(x,'red',attrs=['bold'])
print_blue = lambda x: cprint(x,'blue',attrs=['bold'])

#found this on stack overflow
def split_list(alist, wanted_parts=1):
    length = len(alist)
    return [ alist[i*length // wanted_parts: (i+1)*length // wanted_parts] 
             for i in range(wanted_parts) ]

def skim(lines):
	for line in lines:
		infile = [line]
		outfile = line.replace('.root', '_skim.root').replace('/root','/skim')
		command = "skim/Skim "+infile[0]+" "+outfile
		os.system(command)

os.system("cd current && make && cd ..")
input_file = str(sys.argv[1])
outfile_file = input_file[:-4]
e1d_command = "current/e1d " + input_file + " " + outfile_file.replace('/inputFiles','/outputFiles') + "_test.root"
os.system(e1d_command)

print_red("e1d done")
os.system("cd skim && make && cd ..")
num_cores = multiprocessing.cpu_count()
pool = Pool(processes=num_cores)
lines = [line.rstrip('\n') for line in open(input_file)]

lines_split = split_list(lines, wanted_parts=num_cores)
pool.map(skim, (lines_split))

print_red("Skim done")

input_file = str(sys.argv[1])
outfile_file = input_file[:-4]
e1d_command = "current/e1d " + input_file[:-4] + "_skim.lis" + " " + outfile_file.replace('/inputFiles','/outputFiles') + "_cut.root 2"
os.system(e1d_command)

print_red("Finished\a")
