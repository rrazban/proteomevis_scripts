# [ProteomeVis](https://github.com/rrazban/proteomevis) data curation

The purpose of this repository is to inform the interested ProteomeVis web app user how the data is collected.

Directory and directory trees with an integer as its first character indicate the order in which to run the executable in the corresponding directory. The first executable to be run should be that in `0-identify_structure/0-identify_pdb/`; the second, `0-identify_structure/1-get_pdb/`, and so on.
No numbers within subdirectories indicate that they can be attended to in any order. For example under `1-property_proteomevis/`, `contact_density/`, `ppi_partner/` and `tm_and_sid/` only depend on the chosen PDB structures from executing the scripts in the `0-identify_structure/` subdirectories. Thus, the three subdirectories relevent data files can be generated at any time with respect to eachother.

Some notes:
1. run.py is executed within an organism's directory
- for example, enter `../run.py` in `proteomevis_scripts/0-identify_structure/0-identify_pdb/yeast/`

2. Current files in the repository reflect current data in the ProteomeVis, however, PDB files, structure images and csv files are omitted because of memory. Affected directories which are empty include
- `0-identify_structure/1-get_pdb/'organism'/`
- `0-identify_structure/2-get_pdb_chain/'organism'/`
- `2-proteomevis_backend/*/*/`
