#!/usr/bin/python
import os, sys
from collections import OrderedDict
import subprocess
import datetime
import filecmp

from parse_data import organism, organism_list


today = datetime.date.today()

def tabular_printout(header, output_list, pre_align='', limit = 100):
	if not tabular_printout.counter: print header
	if tabular_printout.counter >= limit: return

	output_list = [str(output) for output in output_list]
	if output_list[0]=="yeast_ecoli": pre_align=''
	print "{0}{1}".format(pre_align, "\t".join(output_list))
	tabular_printout.counter += 1
tabular_printout.counter = 0

def get_full_filename(filename, date_bool):
	if date_bool==True:
		filename+='_{0}.txt'.format(today)
	else: 
		filename+='.txt'
	return filename

def writeout(label_list, d_output, filename='output', date_bool=False):	
	d_output = OrderedDict(sorted(d_output.items()))	
	filename = get_full_filename(filename, date_bool)
	with open(filename, "w") as wfile:
		wfile.write("\t".join(label_list))
		wfile.write("\n")
		for key, value_list in d_output.iteritems():
			if type(value_list)!=list:
				value_list = [value_list]

			wfile.write("{0}\t".format(key))
			if type(value_list[0])==float:
				value_list = [ '%.4f' % value for value in value_list ]
			else: 
				value_list = [ str(value) for value in value_list ]
			wfile.write("\t".join(value_list))
			wfile.write("\n")

def get_next_dirname(CWD, PATH):
	if organism in organism_list:
		CWD = os.path.dirname(CWD)
	current_dirname = os.path.basename(CWD)
	dirname_list = os.listdir(PATH)
	pre_next_dirname = [dirname for dirname in dirname_list if dirname[0].isdigit()]
	next_dirname = [dirname for dirname in pre_next_dirname if int(dirname[0])==int(current_dirname[0])+1]
	if not pre_next_dirname: 
		dirname_list.remove(current_dirname)
		if PATH=='../':
			print "Make sure all scripts in directory are run before moving on"
		print "\nMake sure all other directories {0} are updated before moving on".format(tuple(dirname_list))
	return next_dirname

def print_next_step(PATH='../../'):
	CWD = os.getcwd()
	next_dirname = get_next_dirname(CWD, PATH)
	if len(next_dirname)==1:
		print 'Proceed to ../../{0} to complete update'.format(next_dirname[0])
	elif len(next_dirname)==0:
		next_head_dirname = get_next_dirname(os.path.dirname(CWD), '../{0}'.format(PATH))
		if len(next_head_dirname)==0:
			print 'update complete. Move files to ProteomeVis web app directory'
			print '\tcopy db.sqlite3 into proteomevis/'
			print '\tcopy pdb_image/* and species.*.png into proteomevis/proteomevis/static/img/'
		else:
			print 'Proceed to ../../../{0} to complete update'.format(next_head_dirname[0])
	else: pass
		
def database_update_needed(filename):
        same_file = filecmp.cmp('{0}.txt'.format(filename), '{0}_{1}.txt'.format(filename, today))
	if same_file:
		print "no update needed"
		subprocess.call(['rm', '{0}_{1}.txt'.format(filename, today)])
	else:
		print "updating {0}".format(filename)
		print "'filename' -> 'old_filename'"
		subprocess.call(['mv', '{0}.txt'.format(filename), 'old_{0}.txt'.format(filename)])
		print "'filename-date' -> 'filename'"
		subprocess.call(['mv', '{0}_{1}.txt'.format(filename, today), '{0}.txt'.format(filename)])
		print_next_step()
	return not same_file
