#!/usr/bin/env python
from __future__ import print_function
# Loads ROOT for opening files
from tinydb.middlewares import CachingMiddleware
from tinydb.storages import JSONStorage
from tinydb import TinyDB
import ROOT
from ROOT import gROOT, gBenchmark
import h10
import time
import sys

db = TinyDB("db.json", storage=CachingMiddleware(JSONStorage))

if len(sys.argv) < 2:
    print("Need input root file")
    sys.exit(0)

filename = sys.argv[1]
evnt_id_base = filename.replace("/", ".").split(".")[-2]

data = h10.data(filename)
data.num_entries
start = time.time()
e = 0
for d in data:
    if d.gpart == 0:
        continue
    e += 1
    d.p[0]

    #id_list = []
    # for ids in h10.particle_list[e]:
    #    id_list.append(ids)

    eid = evnt_id_base + str(e)
    #pl = ' '.join(map(str, id_list))
    #pls = ' '.join(map(str, sorted(id_list)))
    #n = len(id_list)

    event = {"evnt_id": eid}
    #         "particle_list": pl,
    #         "particle_list_sorted": pls,
    #         "num_particles": n,
    #         "w": 0,  # h10.W_vec[e],
    #         "q2": 0,  # h10.Q2_vec[e],
    #         "mm": 0}
    db.insert(event)

end = time.time()

print(f'Time {end - start}s => {data.num_entries/(end-start)} Hz')
print(data.num_entries)
