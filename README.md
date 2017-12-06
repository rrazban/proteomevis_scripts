# ProteomeVis data curation

To update the database, run executables (run.py) as 
indicated by the 0, 1, 2, or 3 label at the start 
of the directory name. No number within subdirectory 
indicates that it can be updated in any order. 

run.py is executed within an organism's dir  
ex) enter `../run.py` in `proteomevis_scripts/0-identify_structure/0-identify_pdb/yeast/` 

PDB files, structure images, and the database are 
omitted because of memory. Affected directories:
 - `0-identify_structure/1-get_pdb/'organism'/`
 - `0-identify_structure/2-get_pdb_chain/'organism'/`  
 - `2-proteomevis_backend/*/`
