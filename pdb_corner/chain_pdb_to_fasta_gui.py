# python3
from gooey import *
from Bio.PDB import *
import sys
# input parameters
@Gooey(required_cols=5, program_name='pdb chain to fasta', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser(description="subsets a pdb file by selecting the model and chain and convert the output to fasta format")
    ap.add_argument("-pdb", "--pdb", required=True, widget='FileChooser', help="input pdb file")
    ap.add_argument("-model", "--model", required=True, help="model from pdb file to select(integer)")
    ap.add_argument("-chain", "--chain", required=True, help="chain from pdb file to select")
    ap.add_argument("-id", "--id", required=True, help="pdb id of the protein structure")
    ap.add_argument("-fasta", "--fasta", required=True, widget='FileSaver', help="output fasta file")
    args = vars(ap.parse_args())
    #main
    def seq_from_pdb(structure):
        ppb = PPBuilder()
        for pp in ppb.build_peptides(structure):
            print(">"+args['id']+"_"+args['chain'],pp.get_sequence(), sep="\n")

    parser = PDBParser()
    s = parser.get_structure("name", args['pdb'])
    fill = s[int(args['model'])][args['chain']]
    sys.stdout = open(args['fasta'], 'a')
    seq_from_pdb(fill)
    sys.stdout.close()

if __name__ == '__main__':
    main()
