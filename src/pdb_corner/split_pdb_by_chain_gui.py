# python3
from gooey import *
import Bio
from Bio.PDB import *
import warnings
from Bio import BiopythonWarning
# input parameters
@Gooey(required_cols=1, program_name='split pdb by chain', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
	ap = GooeyParser()
	ap.add_argument("-in", "--input", required=True, widget='FileChooser', help="input pdb file with 1 model")
	args = vars(ap.parse_args())
# main
# ignore warnings
	warnings.simplefilter('ignore', BiopythonWarning)
# retrieve pdb id
	pdbid = args['input'].split('.')[0]
# import pdb
	parser = PDBParser()
	s = parser.get_structure("name", args['input'])
	model = s[0]
# change output directory
	for chain in model:
		io = PDBIO()
		io.set_structure(chain)
		io.save(''.join([pdbid,"_",chain.get_id(),".pdb"]))

if __name__ == '__main__':
    main()