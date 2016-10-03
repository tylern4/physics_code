#!/usr/bin/env python

from tinydb import TinyDB, where
from tinydb.storages import JSONStorage
from tinydb.middlewares import CachingMiddleware
import ROOT

db = TinyDB('/Users/tylern/Desktop/data/data.json',storage=CachingMiddleware(JSONStorage))
file_name = '/Users/tylern/data/inputFiles/v3/skim/r22855_00_skim.root'
file_h10 = ROOT.TFile.Open(file_name, 'read') # a file
num = 0 

for event in file_h10.h10:   # loading the tree
	data = {}
	data['evnt_id'] = file_name.split('/')[-1].split('.')[0] + '_' + str(num)
	num += 1
	data['W'] = event.W
	data['Q2'] = event.Q2
	data['gpart'] = event.gpart
	#data['id'] = [ids for ids in event.id]  
	#data['stat'] = [x for x in str(event.stat)]  
	#data['dc'] = [x for x in str(event.dc)]   
	#data['cc'] = [x for x in str(event.cc)]   
	#data['sc'] = [x for x in str(event.sc)]   
	#data['ec'] = [x for x in str(event.ec)]   
	#data['lec'] = [x for x in str(event.lec)]   
	data['p'] = [x for x in event.p]   
	data['m'] = [x for x in event.m]   
	#data['q'] = [x for x in event.q]   
	data['b'] = [x for x in event.b]   
	data['cx'] = [c for c in event.cx]
	data['cy'] = [c for c in event.cy]
	data['cz'] = [c for c in event.cz]
	data['vx'] = [x for x in event.vx]   
	data['vy'] = [x for x in event.vy]   
	data['vz'] = [x for x in event.vz]
	db.insert(data)
