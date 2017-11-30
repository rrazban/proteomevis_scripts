scripts used to generate the csv files making up the 
ProteomeVis web app (github.com/rrazban/proteomevis/make_database)


order in which executables (run.py) should be run 
is indicated by the 0, 1, or 2 label at the start 
of the directory name. Updates ProteomeVis when new 
PDB structures are found for Uniprot sequences
 - if no number is present, executables (called run.py) can
	be run in any order

update property_proteome/ when proteome is updated

run.py is executed within an organism's dir  
ex) in 0-identify_structure/0-identify_pdb/yeast/, enter ../run.py

PDB files and structure images are omited because of memory
 - 0-identify_structure/1-get_pdb/  
 - 0-identify_structure/2-get_pdb_chain/  
 - 2-make_image/pdb_image/
