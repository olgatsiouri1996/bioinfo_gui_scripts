# python3
from gooey import *
import os
import sys
import argparse
from pyfaidx import Fasta
# input parameters
@Gooey(required_cols=2, program_name= 'extract or remove sequences from fasta', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser(description="use a txt file with fasta headers to extract or remove sequences from fasta file")
    ap.add_argument("-in", "--input", required=True, widget='FileChooser', help="input multi-fasta file")
    ap.add_argument("-ids", "--ids", required=True, widget='FileChooser', help="file with fasta headers to retrieve the output fasta sequences")
    ap.add_argument("-out", "--output", required=False, widget='FileSaver', help="output multi-fasta file")
    ap.add_argument("-dir", "--directory", required=False, type=str, widget='DirChooser',  help="output directory to save the single-fasta files.")
    ap.add_argument("-pro", "--program",type=int, default=1, required=False, help="choose to: 1) extract sequences from a multi-fasta file, 2) extract many single-fasta files, 3) remove sequences from a multi-fasta file, 4) remove sequences and export to many single-fasta files . Defaults to 1")
    args = vars(ap.parse_args())
    # main
    # create function to split the input sequence based on a specific number of characters(60)
    def split_every_60(s): return [str(s)[i:i+60] for i in range(0,len(str(s)),60)]
    # import the txt file with headers you want to extract the sequence from the input fasta
    with open(args['ids'], 'r') as f:
        headers = f.readlines()
    headers = [x.strip() for x in headers]
    # import fasta file
    features = Fasta(args['input'])
    # import program number
    program = args['program']
    # choose program
    match program:
        # extract a multi-fasta file
        case 1:
            # iterate input headers to extract sequences and export as multi-fasta
            sys.stdout = open(args['output'], 'a')
            for header in headers:
                print(''.join([">",features[str(header)].long_name]))
                print('\n'.join(split_every_60(features[str(header)][:].seq)))
            sys.stdout.close()
        case 2:
            # extract many single fasta files
            os.chdir(args['directory'])
            for header in headers:
                sys.stdout = open(''.join([str(header),".fasta"]), 'a')
                print(''.join([">",features[str(header)].long_name]))
                print('\n'.join(split_every_60(features[str(header)][:].seq)))
                sys.stdout.close()
        case 3:
            # remove ids
            keyslist = list(features.keys())
            for header in headers:
                keyslist.remove(header)
            # export to 1 multi-fasta
            sys.stdout = open(args['output'], 'a')
            for key in keyslist:
                print(''.join([">",features[str(key)].long_name]))
                print('\n'.join(split_every_60(features[str(key)][:].seq)))
            sys.stdout.close()
        case 4:
            # remove ids
            keyslist = list(features.keys())
            for header in headers:
                keyslist.remove(header)
            # extract many sigle-fasta files
            os.chdir(args['directory'])
            for key in keyslist:
                sys.stdout = open(''.join([str(key),".fasta"]), 'a')
                print(''.join([">",features[str(key)].long_name]))
                print('\n'.join(split_every_60(features[str(key)][:].seq)))
                sys.stdout.close()

if __name__ == '__main__':
    main()