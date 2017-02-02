#!/usr/local/bin/ipython
from convert import convert, split_list
from StringIO import StringIO
import sys
import os

lines = [line.rstrip('\n') for line in open(str(sys.argv[1]))]

convert(lines)

print("done")