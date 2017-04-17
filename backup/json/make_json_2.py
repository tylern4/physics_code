#!/usr/env python

from tinydb import TinyDB, where
from tinydb.storages import JSONStorage
from tinydb.middlewares import CachingMiddleware
import ROOT

max_value=53629.0
db = TinyDB('test.json',storage=CachingMiddleware(JSONStorage))
num = 0 
file = ROOT.TFile.Open('/mnt/Data/e1d/v2/skim/r23033_08_skim.root', 'read') # a file
for event in file.h10:
	num += 1
	if(num % 540 == 0):
		print (num/max_value)*100

	db.insert({'W': event.W, 'Q2': event.Q2})
