#!/usr/bin/python
import sys, subprocess
import datetime
import filecmp


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
		
def database_update_needed(filename):
        same_file = filecmp.cmp('{0}.txt'.format(filename), '{0}_{1}.txt'.format(filename, today))
	if same_file:
		print "no update needed"
		subprocess.call(['rm', '{0}_{1}.txt'.format(filename, today)])
		return False
	else:
		print "database update needed"
		print "manually change name 'filename-date' -> 'filename'"
		return True
