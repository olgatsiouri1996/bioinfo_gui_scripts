# python3
import sys
from gooey import *
from pyfaidx import Fasta
# input parameters
@Gooey(required_cols=3, program_name='sort multi-fasta', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser()
    ap.add_argument("-in", "--input", required=True, widget='FileChooser', help="input multi-fasta file")
    ap.add_argument("-ids", "--ids", required=True, widget='FileChooser', help="file with sorted fasta identifiers to retrieve the output fasta sequences")
    ap.add_argument("-out", "--output", required=True, widget='FileSaver', help="output multi-fasta file")
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
    # iterate input headers to reorder sequences and export as multi-fasta
    sys.stdout = open(args['output'], 'a')
    for header in headers:
        print(''.join([">",features[str(header)].long_name]))
        print('\n'.join(split_every_60(features[str(header)][:].seq)))
    sys.stdout.close()

if __name__ == '__main__':
    main()