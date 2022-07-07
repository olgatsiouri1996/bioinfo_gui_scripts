# python3
import os
from gooey import *
import Bio
from Bio.PDB import *
import warnings
from Bio import BiopythonWarning
# input parameters
@Gooey(required_cols=3, program_name='subset pdb', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser(description="subsets 1 or many pdb file by selecting the chain and residues from it")
    ap.add_argument("-pdb", "--pdb", required=False, type=str, widget='DirChooser', help=" directory containing the input pdb files")
    ap.add_argument("-in", "--input", required=False, widget='FileChooser', help="input pdb file with 1 model")
    ap.add_argument("-pro", "--program", required=False,default=1, type=int, help="program to choose 1) subset 1 pdb file 2) subset many pdb files. Default is 1")
    ap.add_argument("-chain", "--chain", required=True, help="chain from pdb file to select")
    ap.add_argument("-start", "--start", required=True, type=int, help="amino acid in chain to start writing the pdb file")
    ap.add_argument("-end", "--end", required=True, type=int, help="amino acid in chain to end writing the pdb file")
    args = vars(ap.parse_args())
    #main
    # ignore warnings
    warnings.simplefilter('ignore', BiopythonWarning)
    # choose program
    if args['program'] == 1:
        parser = PDBParser()
        s = parser.get_structure("name", args['input'])
        Bio.PDB.Dice.extract(s, args['chain'], args['start'],args['end'], ''.join([args['input'].split(".")[0],"_",args['chain'],"_",str(args['start']),"_",str(args['end']),".pdb"]))
    else:
        # import each fasta file from the working directory
        for filename in sorted(os.listdir(os.chdir(args['pdb']))):
            if filename.endswith(".pdb"):
                parser = PDBParser()
                s = parser.get_structure("name", filename)
                Bio.PDB.Dice.extract(s, args['chain'], args['start'],args['end'], ''.join([filename.split(".")[0],"_",args['chain'],"_",str(args['start']),"_",str(args['end']),".pdb"]))

if __name__ == '__main__':
    main()


