import os, sys
from collections import OrderedDict
import subprocess
import filecmp

from parse_data import organism, organism_list


def writeout(label_list, d_output, filename='output'):	
	d_output = OrderedDict(sorted(d_output.items()))	
	filename = filename + '.txt' 
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
	if organism in organism_list or organism=='yeast_ecoli':
		CWD = os.path.dirname(CWD)
	current_dirname = os.path.basename(CWD)
	dirname_list = os.listdir(PATH)
	pre_next_dirname = [dirname for dirname in dirname_list if dirname[0].isdigit()]
	next_dirname = [dirname for dirname in pre_next_dirname if int(dirname[0])==int(current_dirname[0])+1]
	if not pre_next_dirname: 
		dirname_list.remove(current_dirname)
		if PATH=='../':
			print "Make sure all scripts in directory are run before moving on"
		print "\nMake sure all other directories {0} are updated before moving on".format(tuple(dirname_list))		#better way: check if files modified after seq2struc.txt using time.ctime(os.path.getmtime(file)))	#also, this way, can catch users starting at 1- or 2- without first doing 0-
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
		
def database_update_needed(filename, bool_print_next_step=True):
        same_file = filecmp.cmp('{0}.txt'.format(filename), 'new_{0}.txt'.format(filename))
	if same_file:
		print "no update needed"
		subprocess.call(['rm', 'new_{0}.txt'.format(filename)])
	else:
		print "updating {0}".format(filename)
		print "'filename' -> 'old_filename'"
		subprocess.call(['mv', '{0}.txt'.format(filename), 'old_{0}.txt'.format(filename)])
		subprocess.call(['mv', 'new_{0}.txt'.format(filename), '{0}.txt'.format(filename)])
		if bool_print_next_step:	#to control TM_SID case	
			print_next_step()

	return not same_file
