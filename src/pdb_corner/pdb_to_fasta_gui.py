# python3
import os
from gooey import *
import warnings
from Bio import BiopythonWarning
from Bio import SeqIO
# input parameters
@Gooey(required_cols=0, program_name='pdb to fasta', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser(description="converts 1 or many pdb files into a fasta file with 1 fasta sequence per chain")
    ap.add_argument("-in", "--input", required=False, widget='FileChooser', help="input pdb file with HEADER section e.g the 1tii structure from Protein Data Bank")
    ap.add_argument("-dir", "--directory", required=False, type=str, widget='DirChooser', help="directory to search for pdb files")
    ap.add_argument("-type", "--type", required=False,default=1, type=int, widget='Dropdown', choices=[1,2],  help="type of input to choose: 1) 1 pdb file, 2) many pdb files.")
    args = vars(ap.parse_args())
    # main
    # ignore warnings
    warnings.simplefilter('ignore', BiopythonWarning)
    # select between importing 1 or many pdb files
    if args['type'] == 1:    
        count = SeqIO.convert(args['input'], "pdb-atom", ''.join([str(args['input']).split('.')[0],'.fasta']), "fasta")
    else:
         # import each pdb file from the working directory
        for filename in sorted(os.listdir(os.chdir(args['directory']))):
            if filename.endswith(".pdb"):
                count = SeqIO.convert(filename, "pdb-atom", ''.join([str(filename).split('.')[0],'.fasta']), "fasta")



if __name__ == '__main__':
    main()