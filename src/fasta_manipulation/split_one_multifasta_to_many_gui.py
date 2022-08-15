# python3
import os
import sys
from gooey import *
from pyfaidx import Fasta
# input parameters
@Gooey(required_cols=2, program_name='split one multi-fasta to many', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser()
    ap.add_argument("-in", "--input", required=True, widget='FileChooser', help="input multi-fasta file")
    ap.add_argument("-num", "--number", required=True, type=int, help="number of fasta records per output fasta file(you can put any number you want as it makes sure the the remainig fasta records will be written to a seperate file as well)")
    args = vars(ap.parse_args())
    # main
    # create function to split the input sequence based on a specific number of characters(60)
    def split_every_60(s): return [str(s)[i:i+60] for i in range(0,len(str(s)),60)]
    # import fasta file
    features = Fasta(args['input'])
    # initial count of files
    count = 0
    # split list
    keyslist = list(features.keys())
    split_lists = [keyslist[x:x+args['number']] for x in range(0, len(keyslist), args['number'])]
    # extract many sigle-fasta files
    for lis in split_lists:
        count = count + 1
        sys.stdout = open(''.join([str(args['input']).split('.fa')[0],"_","part",str(count),".fasta"]), 'a')
        for key in lis:
            print(''.join([">",features[str(key)].long_name]).replace('\r',''))
            print('\n'.join(split_every_60(features[str(key)][:].seq)))
        sys.stdout.close()

if __name__ == '__main__':
    main()