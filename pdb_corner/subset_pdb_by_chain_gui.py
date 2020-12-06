# python3
from gooey import *
from Bio.PDB import *
# input parameters
@Gooey(required_cols=3, program_name='subset pdb by chain', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser(description="subsets a pdb file by selecting the model and chain from it")
    ap.add_argument("-in", "--input", required=True, widget='FileChooser', help="input pdb file")
    ap.add_argument("-model", "--model", required=False, default= 0, help="model from pdb file to select(integer, default=0)")
    ap.add_argument("-chain", "--chain", required=True, help="chain from pdb file to select")
    ap.add_argument("-out", "--output", required=True, widget='FileSaver', help="output pdb file")
    args = vars(ap.parse_args())
    # main
    parser = PDBParser()
    s = parser.get_structure("name", args['input'])
    fill = s[int(args['model'])][args['chain']]
    io = PDBIO()
    io.set_structure(fill)
    io.save(args['output'])

if __name__ == '__main__':
    main()