# python3
import os
from gooey import *
import Bio
from Bio.PDB import *
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
# setup empty lists
    ids = []
    chains = []
    star_pos = []
    end_pos = []    
# inport txt file and convert each column to list
    with open(args['txt'], 'r') as f:
            for line in f:
                # convert each column to list
                ids.append(line.split()[0])
                chains.append(line.split()[1])
                star_pos.append(line.split()[2])
                end_pos.append(line.split()[3])
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
