# python3
from gooey import *
from Bio.PDB import *
import pandas as pd
import os
# input parameters
@Gooey(required_cols=3, program_name='dssp statistics by chain', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser(description="select model and chain from pdb and output secondary strucure statistics using DSSP")
    ap.add_argument("-in", "--input_pdb", required=True, widget='FileChooser', help="input pdb file")
    ap.add_argument("-model", "--pdb_model", required=False, default= 0, help="model from pdb file to select(integer, default=0)")
    ap.add_argument("-chain", "--pdb_chain", required=True, help="chain from pdb file to select")
    ap.add_argument("-out", "--output", required=True, widget='FileSaver', help="output txt file(every name except: out.pdb)")
    args = vars(ap.parse_args())
# main
# select model and chain
    parser = PDBParser()
    s = parser.get_structure("name", args['input_pdb'])
    fill = s[int(args['pdb_model'])][args['pdb_chain']]
    io = PDBIO()
    io.set_structure(fill)
    io.save("out.pdb")
# import trimed pdb file and run dssp 
    parser = PDBParser()
    s = parser.get_structure("name", "out.pdb")
    fill = s[0]
    dssp = DSSP(fill,"out.pdb", dssp='mkdssp')
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
# remove intermediate file
    os.system("rm *out.pdb")

if __name__ == '__main__':
    main()