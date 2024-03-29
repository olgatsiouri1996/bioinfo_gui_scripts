# python3
from gooey import *
import sys
from pyfaidx import Fasta
import  pandas as pd
# input parameters
@Gooey(required_cols=3, program_name='trim protein by coords', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser(description='trim a multi-fasta file with proteins based on a txt file with start and end coordinates')
    ap.add_argument("-in", "--input", required=True, widget="FileChooser", help="input fasta file")
    ap.add_argument("-coords", "--coordinates", required=True, widget="FileChooser", help="input 3-column tab-seperated txt file with id, start and end positions respectively in each row")
    ap.add_argument("-out", "--output", required=True, widget="FileSaver", help="output multi-fasta file")
    args = vars(ap.parse_args())
    # main
    # create function to split the input sequence based on a specific number of characters(60)
    def split_every_60(s): return [str(s)[i:i+60] for i in range(0,len(str(s)),60)]
    # inport txt file and convert each column to list
    df_txt = pd.read_csv(args['coordinates'], header=None, sep="\t")
    ids = df_txt.iloc[:,0].values.tolist()
    seq_start = df_txt.iloc[:,1].values.tolist()
    seq_start[:] = [i - 1 for i in seq_start]
    seq_end = df_txt.iloc[:,2].values.tolist()
    # import fasta file
    features = Fasta(args['input'])
    # iterate all below lists in pairs
    sys.stdout = open(args['output'], 'a')
    for (a, b, c) in zip(ids, seq_start, seq_end):
        print(''.join([">",str(a),"_",str(b),"_",str(c)]))
        print('\n'.join(split_every_60(features[str(a)][int(b):int(c)].seq)))
    sys.stdout.close()

if __name__ == '__main__':
    main()
