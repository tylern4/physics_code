#!/usr/bin/env python
from __future__ import print_function
# Loads ROOT for opening files
import ROOT
from ROOT import gROOT, gBenchmark

import sys
import cppyy
cppyy.load_reflection_info("H10.so")

from pymongo import MongoClient
client = MongoClient('mongodb://172.21.139.37:27017')
db = client['e1d']
events = db['events']

if len(sys.argv) < 2:
    print("Need input root file")
    sys.exit(0)

filename = sys.argv[1]
evnt_id_base = filename.replace("/", ".").split(".")[-2]

gROOT.SetBatch(True)
chain = ROOT.TChain('h10')
h10 = cppyy.gbl.H10()

num_files = chain.Add(filename)
h10_data = h10.dataHandeler(chain)

num = len(h10_data.W_vec)

for e in xrange(0, num):
    id_list = []
    for ids in h10_data.particle_list[e]:
        id_list.append(ids)

    eid = evnt_id_base + "_" + str(e)
    pl = ' '.join(map(str, id_list))
    pls = ' '.join(map(str, sorted(id_list)))
    n = len(id_list)

    event = {"evnt_id": eid,
             "particle_list": pl,
             "particle_list_sorted": pls,
             "num_particles": n,
             "w": h10_data.W_vec[e],
             "q2": h10_data.Q2_vec[e],
             "mm": 0}
    events.insert_one(event)


print(evnt_id_base)
