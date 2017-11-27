scripts used to generate the csv files making up the 
ProteomeVis web app (github.com/rrazban/proteomevis/make_database)

order in which executables (run.py) should be run to 
is indicated by the 0, 1, or 2 label at the start of 
the directory name

if no number present, executables can be run in any order

run.py is run within an organism's dir
ex) in 0-identify_structure/0-identify_pdb/yeast/, enter ../run.py

PDB files and structure images are omited because of memory
considerations
PDB files:	0-identify_structure/1-get_pdb/*/
		0-identify_structure/2-get_pdb_chain/*/
PDB image:	2-make_image/pdb_image/*
