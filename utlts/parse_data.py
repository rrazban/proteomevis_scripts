import os, sys	#call module 'organism_setup?'


organism_list = ['yeast','ecoli']#, 'protherm']#,'yeast_ecoli' append when wanted	#need to figure out a better way

cwd = os.getcwd()
organism = os.path.basename(cwd)

def initialize_dict(dtype):	#not sure if this belongs in this module
	d = {}

	for organism in organism_list:
		if dtype=="list":
			d[organism] = []
		elif dtype=="dict":
			d[organism] = {} 
		else:
			print "unavailable type for dict formation"
			sys.exit()
	return d

def taxid():
	d = {}
	d['yeast'] = 559292
	d['ecoli'] = 83333
	return d

def int2organism():
	d={}
	d[0] = 'yeast'
	d[1] = 'ecoli'
	return d
