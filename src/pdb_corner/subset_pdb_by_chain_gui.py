# python3
from gooey import *
import Bio
from Bio.PDB import *
import warnings
from Bio import BiopythonWarning
# input parameters
@Gooey(required_cols=3, program_name='subset pdb by chain', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
	ap = GooeyParser()
	ap.add_argument("-in", "--input", required=True, widget='FileChooser', help="input pdb file with 1 model")
	ap.add_argument("-chain", "--chain", required=True, help="chain from pdb file to select")
	ap.add_argument("-out", "--output", required=True, widget='FileSaver', help="output pdb file")
	args = vars(ap.parse_args())
# main
	warnings.simplefilter('ignore', BiopythonWarning)
	parser = PDBParser()
	s = parser.get_structure("name", args['input'])
	model = s[0]
	chain = model[args['chain']]
	io = PDBIO()
	io.set_structure(chain)
	io.save(args['output'])

if __name__ == '__main__':
    main()