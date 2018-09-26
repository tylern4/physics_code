#!/usr/bin/env python
import glob
import sys
import time
from multiprocessing import Pool

import skim


class colors:
    HEADER = '\033[95m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    WARNING = '\033[93m'
    RED = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


USAGE = """
{} /path/to/input/files/ output.root
"""

if len(sys.argv) == 1:
    print(colors.WARNING + "No inputs given: " + colors.BLUE + USAGE.format(sys.argv[0]) + colors.ENDC)
    sys.exit(2)
elif len(sys.argv) == 2:
    input_files = sys.argv[1]
    output_file = "out.root"
elif len(sys.argv) == 3:
    input_files = sys.argv[1]
    output_file = sys.argv[2]

files = glob.glob(input_files + "*.root")

if sys.version_info[0] > 3:
    map(str.encode, files)
    output_file = str.encode(output_file)


start = time.time()
for f in files:
    skim.py_skim([f], f + "_skim.root")

end = time.time()
print(colors.GREEN + "{} Sec".format(end - start) + colors.ENDC)
