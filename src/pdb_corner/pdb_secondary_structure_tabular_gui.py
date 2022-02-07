#python3
from gooey import *
from Bio.PDB import *
import pandas as pd
# input parameters
@Gooey(required_cols=2, program_name='pdb secondary stracture statistics', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser(description="output a 3 column tabular file with residue number, name and secondary strucure using DSSP")
    ap.add_argument("-pdb", "--input", required=True, widget='FileChooser', help="input pdb file")
    ap.add_argument("-model", "--model", required=False, default= 0, help="model from pdb file to select(integer, default=0)")
    ap.add_argument("-out", "--output", required=True, widget='FileSaver', help="output txt file")
    args = vars(ap.parse_args())
    # main
    parser = PDBParser()
    s = parser.get_structure("name", args['input'])
    fill = s[int(args['model'])]
    dssp = DSSP(fill, args['input'], dssp='mkdssp')
    df = pd.DataFrame(dssp)
    df = df.loc[:, [0,1,2]]
    # export
    with open(args['output'], 'a') as f:
        f.write(
            df.to_csv(header = ['residue_number','residue_name', 'residue_structure'], index = False, sep= "\t", doublequote= False, line_terminator= '\n')
        )

if __name__ == '__main__':
    main()