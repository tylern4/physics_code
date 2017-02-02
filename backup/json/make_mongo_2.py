#!/usr/bin/env python
from pymongo import MongoClient
import ROOT

client = MongoClient('mongodb://172.21.139.25:27017/')
db = client.physics
collection = db.v3
phyDB = db.posts
num = 0

file = ROOT.TFile.Open('/Users/tylern/data/inputFiles/v3/skim/r22855_01_skim.root', 'read') # a file
for event in file.h10:   # loading the tree
	phyDB.insert_one({'W': event.W, 'Q2': event.Q2})
	num += 1
	print num
