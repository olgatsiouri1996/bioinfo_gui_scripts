# python3
import os
from gooey import *
import Bio
from Bio.PDB import *
import warnings
from Bio import BiopythonWarning
# input parameters
@Gooey(required_cols=1, program_name='subset pdb by chain', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
	ap = GooeyParser()
	ap.add_argument("-in", "--input", required=False, widget='FileChooser', help="input pdb file with 1 model")
	ap.add_argument("-pro", "--program", required=False,default=1, type=int, help="program to choose 1) subset 1 pdb file 2) subset many pdb files. Default is 1")
	ap.add_argument("-pdb", "--pdb", required=False, type=str, widget='DirChooser', help=" directory containing the input pdb files")
	ap.add_argument("-chain", "--chain", required=True, help="chain from pdb file to select")
	args = vars(ap.parse_args())
	# main
	# ignore warnings
	warnings.simplefilter('ignore', BiopythonWarning)
	# choose program
	if args['program'] == 1:
	    parser = PDBParser()
	    s = parser.get_structure("name", args['input'])
	    model = s[0]
	    chain = model[args['chain']]
	    io = PDBIO()
	    io.set_structure(chain)
	    io.save(''.join([args['input'].split(".")[0],"_",args['chain'],".pdb"]))
	else:
	    # import each fasta file from the working directory
	    for filename in sorted(os.listdir(os.chdir(args['pdb']))):
	        if filename.endswith(".pdb"):
	    	    parser = PDBParser()
	    	    s = parser.get_structure("name", filename)
	    	    model = s[0]
	    	    chain = model[args['chain']]
	    	    io = PDBIO()
	    	    io.set_structure(chain)
	    	    io.save(''.join([filename.split(".")[0],"_",args['chain'],".pdb"]))

if __name__ == '__main__':
	main()
