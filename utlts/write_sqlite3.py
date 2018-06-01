import sqlite3

class SQLite3():
	def __init__(self, table_name, columns, line_list):
		self.table_name = table_name
		self.columns = columns
		self.line_list = line_list

	def connect(self):
		self.conn = sqlite3.connect('db.sqlite3')
		self.c = self.conn.cursor()
		self.c.execute('DROP TABLE IF EXISTS {0}'.format(self.table_name))
		self.c.execute('CREATE TABLE {0}({1})'.format(self.table_name, ','.join(self.columns)))

	def write(self):
		for line in self.line_list:
			self.c.execute("INSERT INTO {0} VALUES {1}".format(self.table_name, tuple(line)))

	def handle_null(self):
	        for column in self.columns: #unable to pass NULL through python list
	                self.c.execute("UPDATE {0} SET {1}=null where {1}=''".format(self.table_name, column))

	def close(self):
		self.conn.commit()
		self.conn.close()

	def run(self):
		self.connect()
		self.write()
		self.handle_null()
		self.close()
