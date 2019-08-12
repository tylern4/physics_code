#!/usr/bin/env python
from __future__ import print_function, unicode_literals, absolute_import
# Loads ROOT for opening files
from pymongo import MongoClient
import sys
import time
import h10
import ROOT
from ROOT import gROOT, gBenchmark

from halo import Halo
#{'text': 'Startup', 'spinner': 'dots'}


client = MongoClient('mongodb://172.21.139.204:27017')
print(client.list_database_names())
db = client['e1d']
events = db['events']

if len(sys.argv) < 2:
    print("Need input root file")
    sys.exit(0)

filename = sys.argv[1]
evnt_id_base = filename.replace("/", ".").split(".")[-2]


spinner = Halo()
spinner.succeed('Starting DataHandeler')
spinner.start('Starting DataHandeler')

spinner.start('Filling mongodb')

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
    db.events.insert_one(event)

end = time.time()

print(f'Time {end - start}s => {data.num_entries/(end-start)} Hz')
print(data.num_entries)

spinner.stop()
