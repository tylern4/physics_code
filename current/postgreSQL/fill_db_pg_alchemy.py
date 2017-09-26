#!/usr/bin/env python
from __future__ import print_function
# Loads ROOT for opening files
import ROOT
from ROOT import gROOT, gBenchmark

import sys
import cppyy
cppyy.load_reflection_info("H10.so")


from sqlalchemy import Table, Column, Integer, String, ForeignKey


def connect(db, host='localhost', port=5432):
    import sqlalchemy
    '''Returns a connection and a metadata object'''
    # We connect with the help of the PostgreSQL URL
    # postgresql://federer:grandestslam@localhost:5432/tennis
    url = 'postgresql://{}@{}:{}/{}'
    url = url.format("postgres", host, port, db)

    # The return value of create_engine() is our connection object
    con = sqlalchemy.create_engine(url, client_encoding='utf8')

    # We then bind the connection to MetaData()
    meta = sqlalchemy.MetaData(bind=con, reflect=True)

    return con, meta


con, meta = connect('e1d', host='172.21.139.37')
events = meta.tables['events']

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

    clause = events.insert().values(evnt_id=eid,
                                    particle_list=pl,
                                    particle_list_sorted=pls,
                                    num_particles=n,
                                    w=h10_data.W_vec[e],
                                    q2=h10_data.Q2_vec[e],
                                    mm=0)
    result = con.execute(clause)


print(evnt_id_base)
