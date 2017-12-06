#!/usr/bin/python

help_msg = 'generate all possible UniProt sequences with matching PDB structures using the SIFTS resource'	#selection criterion 1
#UniProt is updated on a rolling basis
#PDB and SIFTS is updated every Tuesday at 7ET (0 UTC on Wednesday)
#code adapted from PDBe example: https://github.com/PDBeurope/PDBe_Programming

import os, sys
import urllib2
import json

CWD = os.getcwd()
UTLTS_DIR = CWD[:CWD.index('proteomevis_scripts')]+'/proteomevis_scripts/utlts'
sys.path.append(UTLTS_DIR)
from parse_user_input import help_message
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

def get_all_info(total_batch, uniprot_list):
	full_url = "https://www.ebi.ac.uk/pdbe/api/mappings/best_structures"
	increment = 1000	#1000 limit set by PDBe	#in current setup, if no structures in range, program fails
	d = {}

	for i in range(total_batch/increment):
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

def get_best_pdb_chain(d_uniprot_info):
	g_coverage=0	#g=degeneracy
	g_resolution=0
	d = {}

	for uniprot, info_list in d_uniprot_info.iteritems():
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
		d[uniprot.encode('utf-8')] = pdb_chain
	return d, g_coverage, g_resolution 

def print_verbose(proteome_size, pre_seq2struc_size, g_coverage, g_resolution):
	print "Proteome size: {0}".format(proteome_size)
	print "# proteins in proteome with structure: {0} ({1:.2f})".format(pre_seq2struc_size, pre_seq2struc_size/float(proteome_size))
	print "\tdegeneracy of those with the same max coverage: {0} ({1:.2f})".format(g_coverage, g_coverage/float(pre_seq2struc_size))
	print "\t\tdegeneracy of those also with the same max resolution: {0} ({1:.2f})".format(g_resolution, g_resolution/float(pre_seq2struc_size))


def prepare_writeout(d_uniprot_pdb, d_proteome):
	d_output = {}
	for uniprot, pdb in d_uniprot_pdb.iteritems():
		d_output[uniprot] = [pdb, d_proteome[uniprot]]
	return d_output
	

if __name__ == '__main__':
	args = help_message(help_msg, bool_add_verbose=True)
	d_proteome = read_in('Entry', 'Gene names  (ordered locus )', filename='proteome')	
	
	proteome_size = len(d_proteome)
	d_uniprot_info = get_all_info(proteome_size, d_proteome.keys())
	d_uniprot_pdb, g_coverage, g_resolution = get_best_pdb_chain(d_uniprot_info)
	pre_seq2struc_size = len(d_uniprot_pdb)

	if args.verbose:
		print_verbose(proteome_size, pre_seq2struc_size, g_coverage, g_resolution)
	d_output = prepare_writeout(d_uniprot_pdb, d_proteome)
	filename='pre_seq2struc'
	writeout(['uniprot', 'pdb', 'oln'], d_output, filename=filename, date_bool=True)
	database_update_needed(filename=filename)
