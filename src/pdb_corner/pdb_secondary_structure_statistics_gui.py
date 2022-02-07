#python3
from gooey import *
from Bio.PDB import *
import pandas as pd
# input parameters
@Gooey(required_cols=2, program_name='pdb secondary stracture statistics', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser(description="output a 2 column tabular file with the secondary strucure and its percentange using DSSP")
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
    df = df.loc[:, 2]
    struct_list = df.values.tolist()
    df1 = pd.DataFrame()
    df1['struct_list'] = struct_list
    df1 = df1['struct_list'].value_counts()
    df1 = round((df1 / df1.sum(axis=0)) * 100, 2)
    # export
    with open(args['output'], 'a') as f:
        f.write(
            df1.to_csv(header = False, index = True, sep= "\t", doublequote= False, line_terminator= '\n')
        )

if __name__ == '__main__':
    main()