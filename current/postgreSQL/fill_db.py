#!/usr/bin/env python
from __future__ import print_function
# Loads ROOT for opening files
import ROOT
from ROOT import gROOT, gBenchmark

import sys
import cppyy
cppyy.load_reflection_info("H10.so")

import psycopg2

sql = """INSERT INTO events (evnt_id,particle_list,particle_list_sorted,num_particles,w,q2,mm)
        VALUES (%s,%s,%s,%s,%s,%s,%s);"""
connection = psycopg2.connect(
    host="172.21.139.37", database="e1d", user="postgres")


class Record():

    def __init__(self, connection):
        self.connection = connection
        self.cursor = connection.cursor()
        self.evnt_id = None
        self.particle_list = None
        self.particle_list_sorted = None
        self.num_particles = 0
        self.w = 0
        self.q2 = 0
        self.mm = 0
        self.sql = None

    def __repr__(self):
        return "Record()"

    def __str__(self):
        if self.sql is None:
            string = """(%s,%s,%s,%s,%s,%s,%s)""" % (self.evnt_id, self.particle_list,
                                                     self.particle_list_sorted,
                                                     self.num_particles, self.w, self.q2, self.mm)
        else:
            string = self.sql % (self.evnt_id, self.particle_list,
                                 self.particle_list_sorted,
                                 self.num_particles, self.w, self.q2, self.mm)
        return string

    def set_sql(self, sql):
        self.sql = sql

    def fill_record(self, **kwargs):
        self.evnt_id = kwargs.get('evnt_id', self.evnt_id)
        self.particle_list = kwargs.get('particle_list', self.particle_list)
        self.num_particles = kwargs.get('num_particles', self.num_particles)
        self.w = kwargs.get('w', self.w)
        self.q2 = kwargs.get('q2', self.q2)
        self.mm = kwargs.get('mm', self.mm)

    def _check(self):
        if not isinstance(self.particle_list, basestring):
            self.num_particles = len(self.particle_list)
            self.particle_list_sorted = ' '.join(
                map(str, sorted(self.particle_list)))
            self.particle_list = ' '.join(map(str, self.particle_list))

        if self.num_particles == 0:
            return False
        elif self.sql is None:
            print("Remeber to set sql statement")
            return False
        elif self.particle_list is None:
            return False
        else:
            return True

    def _reset(self):
        self.evnt_id = None
        self.particle_list = None
        self.particle_list_sorted = None
        self.num_particles = 0
        self.w = 0
        self.q2 = 0
        self.mm = 0

    def insert_record(self):
        if self._check():
            self.cursor.execute(self.sql, (self.evnt_id, self.particle_list,
                                           self.particle_list_sorted,
                                           self.num_particles, self.w, self.q2, self.mm))
            self.connection.commit()
        self._reset()

# cur.executemany(sql, seq_of_parameters)
# cur.execute(sql,(2,3,'11 221','11 221','11 221',1,0,0,0))
# conn.commit()
# TODO:
# Open root file
# loop over file
# make evnt_id from file name and event number
# Get list of for particle_list into string (join with ' '.join(...))
# Sort the list and put in particle_list_sorted
# For now use particle_list as my_particle_list
# Calcualte W and Q2
# For not put in MM = 0
if len(sys.argv) < 2:
    print("Need input root file")
    sys.exit(0)

filename = sys.argv[1]
evnt_id_base = filename.replace("/", ".").split(".")[-2]


rec = Record(connection)
rec.set_sql(sql)

gROOT.SetBatch(True)
chain = ROOT.TChain('h10')
h10 = cppyy.gbl.H10()

num_files = chain.Add(filename)
h10_data = h10.dataHandeler(chain)

num = len(h10_data.W_vec)

for e in xrange(0, num):
    id_list = []
    rec.fill_record(evnt_id=evnt_id_base + "_" + str(e))
    for ids in h10_data.particle_list[e]:
        id_list.append(ids)

    rec.fill_record(w=h10_data.W_vec[e])
    rec.fill_record(q2=h10_data.Q2_vec[e])
    rec.fill_record(particle_list=id_list)
    rec.insert_record()

"""

CREATE TABLE EVENTS(
EVNT_ID TEXT NOT NULL,
NUM_PARTICLES INT,
PARTICLE_LIST TEXT,
PARTICLE_LIST_SORTED TSVECTOR,
W REAL,
Q2 REAL,
MM REAL);

"""
