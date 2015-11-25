import os

#found this on stack overflow
def split_list(alist, wanted_parts=1):
    length = len(alist)
    return [ alist[i*length // wanted_parts: (i+1)*length // wanted_parts] 
             for i in range(wanted_parts) ]

def skim(lines):
	for line in lines:
		infile = [line]
		outfile = line.replace('.root', '_skim.root').replace('/root','/skim')
		command = "./Skim "+infile[0]+" "+outfile
		os.system(command)
