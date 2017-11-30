#!/usr/bin/python

#selection criterion 1: generate all possible UniProt sequence to PDB matching using SIFTS
#UniProt is updated on a rolling basis
#PDB and SIFTS is updated every Tuesday at 7ET (0 UTC on Wednesday)
#code adapted from PDBe example: https://github.com/PDBeurope/PDBe_Programming

import sys, collections
import urllib2
import json

sys.path.append('../../../utlts/')
from read_in_file import read_in
from output import writeout, database_update_needed


def make_request(url, data):
    request = urllib2.Request(url)
    try:
        url_file = urllib2.urlopen(request, data)
    except urllib2.HTTPError as e:
        if e.code == 404:
            print("[NOTFOUND %d] %s" % (e.code, url))
        else:
            print("[ERROR %d] %s" % (e.code, url))
        return None
    return url_file.read().decode()

def post_request(data, full_url):
	if isinstance(data, (list, tuple)):
		data = ",".join(data)
	return make_request(full_url, data.encode())

def get_all_info(total, uniprot_list):
	full_url = "https://www.ebi.ac.uk/pdbe/api/mappings/best_structures"
	increment = 1000	#1000 limit set by PDBe	#in current setup, if no structures in range, program fails
	d = {}
	for i in range(total/increment):
		output = post_request(uniprot_list[i*increment:(i+1)*increment], full_url)
		d_output = json.loads(output)
		d.update(d_output)
	return d

def get_resolution(new_info_list):
	data = []
	for info in new_info_list:
		data.append(info[u'resolution'])
	return data

def get_coverage(info_list):
	data = []
	for info in info_list:
		coverage = info[u'coverage']	#make sure matches with gene
		data.append(coverage)
	return data

def update_info_list(info_list, index_list):
	data = []
	for index in index_list:
		data.append(info_list[index])
	return data

def get_best_pdb_chain(d_info, d_uni):
	d = {}
	g_coverage=0	#g=degeneracy
	g_resolution=0
	for uniprot, info_list in d_info.iteritems():
		if uniprot=='P0CX34':
			print info_list
		coverage_list = get_coverage(info_list)
		max_coverage_list = [i for i, coverage in enumerate(coverage_list) if coverage == max(coverage_list)]	#account for ties
		if len(max_coverage_list)>1:
			g_coverage+=1
			new_info_list = update_info_list(info_list, max_coverage_list)
			resolution_list = get_resolution(new_info_list)
			min_resolution_list = [i for i, resolution in enumerate(resolution_list) if resolution == min(resolution_list)]
			if len(min_resolution_list) > 1:
				g_resolution+=1
			info = new_info_list[min_resolution_list[0]]	#arbitrarily choose first item in list
		else:
			info = info_list[max_coverage_list[0]]
		pdb_basename = info[u'pdb_id']
		pdb_chain_id = info[u'chain_id']
		pdb_chain = "{0}.{1}".format(pdb_basename, pdb_chain_id)
		d[pdb_chain] = uniprot.encode('utf-8')
	return d, g_coverage, g_resolution 

def add_ensembl(d):
	d_ref = read_in('Entry', 'Gene names  (ordered locus )', filename='proteome')		
	new_d = {}
	for pdb, uniprot in d.iteritems():
		new_d[uniprot] = [pdb, d_ref[uniprot]]
	return new_d
	

if __name__ == '__main__':
	d_uni = read_in('Entry', 'Length', filename='proteome')	
	
	total = len(d_uni)
	print "Proteome size: {0}".format(total)

	uniprot_list = d_uni.keys()	
	d_info = get_all_info(total, uniprot_list)
	d_uni_pdb, g_coverage, g_resolution = get_best_pdb_chain(d_info, d_uni)
	total_struc = len(d_uni_pdb)

	print "# proteins in proteome with structure: {0} ({1:.2f})".format(total_struc, total_struc/float(total))
	print "Degeneracy of those proteins with structures with the same max coverage: {0} ({1:.2f})".format(g_coverage, g_coverage/float(total_struc))
	print "Degeneracy of those proteins with structures with the same max resolution: {0} ({1:.2f})".format(g_resolution, g_resolution/float(total_struc))

	d_uni_list = add_ensembl(d_uni_pdb)
	filename='pre_seq2struc'
	writeout(['uniprot', 'pdb', 'oln'], collections.OrderedDict(sorted(d_uni_list.items())), filename=filename, date_bool=True)

	database_update_needed(filename=filename)
