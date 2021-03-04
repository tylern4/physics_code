#!/usr/bin/env python
from __future__ import print_function
import lxml.html as PARSER
import sys

data = open(sys.argv[1]).read()
root = PARSER.fromstring(data)

trs = []
for tr in root.getiterator():
    if tr.tag == "tr":
        trs.append(tr.text_content())


out = [x.replace("\n", ",") for x in trs]
out = [x.strip('\n').strip(' ') for x in out]
out = [x.replace(", ", ",") for x in out]
out = [x.replace(" ,", ",") for x in out]
out = [x for x in out if ",,," not in x]
out = [x.split(",") for x in out]
out = [','.join(x[0:9]) for x in out]


output = sys.argv[2] if len(sys.argv) > 2 else "out.csv"

with open(output, "w") as f:
    for c in out[1:]:
        print(c)
        print(c, file=f)
