# python3
from gooey import *
import os
import sys
import argparse
from pyfaidx import Fasta
# input parameters
@Gooey(required_cols=3, program_name= 'extract one sequences from fasta', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser()
    ap.add_argument("-in", "--input", required=True, widget='FileChooser', help="input multi-fasta file")
    ap.add_argument("-id", "--identifier", required=True, type=str, help="fasta identifier to retrieve the fasta sequence for")
    ap.add_argument("-dir", "--directory", required=False, type=str, widget='DirChooser',  help="output directory to save the single-fasta file.")
    args = vars(ap.parse_args())
    # main
    # create function to split the input sequence based on a specific number of characters(60)
    def split_every_60(s): return [str(s)[i:i+60] for i in range(0,len(str(s)),60)]
    # import multi-fasta file
    features = Fasta(args['input'])
    # extract 1 fasta record based on a fasta identifier
    os.chdir(args['directory'])
    sys.stdout = open(''.join([args['identifier'],".fasta"]), 'a')
    print(''.join([">",features[args['identifier']].long_name]))
    print('\n'.join(split_every_60(features[args['identifier']][:].seq)))
    sys.stdout.close()

if __name__ == '__main__':
    main()
