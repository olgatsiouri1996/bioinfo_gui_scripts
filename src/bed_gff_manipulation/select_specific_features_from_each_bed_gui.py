# python3
import os
from gooey import *
import pandas as pd
import warnings
# input parameters
@Gooey(required_cols=2, program_name= 'select specific features from each bed file', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser()
    ap.add_argument("-in", "--input", required=True, widget='FileChooser', help="input txt file with filenames(with the \'.bed\' removed)")
    ap.add_argument("-fea", "--feature", required=False, default='gene',  type=str, help="specify the pattern to print the lines that have it")
    ap.add_argument("-col", "--column", required=False, default=8, type=int, help="specify the column number to search for the pattern")
    ap.add_argument("-dir", "--directory", required=True, type=str, widget='DirChooser',  help="directory with bed files")
    args = vars(ap.parse_args())
    # main
    # ignore warnings
    warnings.filterwarnings('ignore')
    # fix index for column
    feature_col = args['column'] - 1
    # import the txt file with filenames
    with open(args['input'], 'r') as f:
        orgs = f.readlines()
    orgs = [x.strip() for x in orgs]
    # go to the directory the .bed files are located
    os.chdir(args['directory'])
    # iterate for each file in the current directory
    for org in orgs:
        # import bed with no headers specified
        df = pd.read_csv(''.join([org,".bed"]), sep= "\t", header=None)
        # select the rows containing the feature
        bool2 = df.iloc[:, feature_col].str.contains(args['feature']) 
        df = df[bool2]
        # export
        with open(''.join([org,"_",args['feature'],".bed"]), 'a') as f:
            f.write(
                df.to_csv(header = False, index = False, sep= "\t", doublequote= False, line_terminator= '\n')
            )

        del df; del bool2

if __name__ == '__main__':
    main()