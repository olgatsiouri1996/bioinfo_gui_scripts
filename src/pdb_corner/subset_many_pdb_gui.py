# python3
import os
from gooey import *
import Bio
from Bio.PDB import *
import  pandas as pd
import warnings
from Bio import BiopythonWarning
# input parameters
@Gooey(required_cols=3, program_name='subset many pdb files', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
	ap = GooeyParser()
	ap.add_argument("-pdb", "--pdb", required=True, type=str, widget='DirChooser', help=" directory containing the input pdb files")
	ap.add_argument("-txt", "--txt", required=True, widget='FileChooser', help="4 column tabular txt file with pdb filename, chain, start and end position columns respectively")
	ap.add_argument("-out", "--out", required=True, type=str, widget='DirChooser', help=" directory to save the output pdb files")
	args = vars(ap.parse_args())
# main
# ignore warnings
	warnings.simplefilter('ignore', BiopythonWarning)
# inport txt file and convert each column to list
	df_txt = pd.read_csv(args['txt'], header=None, sep="\t")
	ids = df_txt.iloc[:,0].values.tolist()
	chains = df_txt.iloc[:,1].values.tolist()
	star_pos = df_txt.iloc[:,2].values.tolist()
	end_pos = df_txt.iloc[:,3].values.tolist()	
# iterate all above lists and subset each pdb file
	for (a, b, c, d) in zip(ids,chains,star_pos,end_pos):
# set working directory
		os.chdir(args['pdb'])
# import pdb
		parser = PDBParser()
		s = parser.get_structure("name", ''.join([str(a),".pdb"]))
# set working directory to export
		os.chdir(args['out'])
# supbset and export to pdb
		Bio.PDB.Dice.extract(s, str(b), int(c), int(d), ''.join([str(a),"_",str(b),"_",str(c),"_",str(d),".pdb"]))

if __name__ == '__main__':
    main()
