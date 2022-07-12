# python3
import sys
from gooey import *
from pyfaidx import Fasta
# imput parameters
@Gooey(required_cols=3, program_name='intersect fasta', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser()
    ap.add_argument("-fa1", "--fasta1", required=True, widget='FileChooser', help="input multi fasta file")
    ap.add_argument("-fa2", "--fasta2", required=True, widget='FileChooser', help="input multi fasta file")
    ap.add_argument("-out", "--output", required=True, widget='FileSaver',  help="output multi fasta file")
    ap.add_argument("-pro", "--program", type=int, default=1, required=False, help="program to choose: 1. export the sequences of fasta1 that have the same identifier as fasta2, 2. export the sequences of fasta1 that doesn't have the same identifier as fasta2. Defaults to 1")
    args = vars(ap.parse_args())
    # main
    # create function to split the input sequence based on a specific number of characters(60)
    def split_every_60(s): return [str(s)[i:i+60] for i in range(0,len(str(s)),60)]
    # create fasta index
    features1 = Fasta(args['fasta1'])
    features2 = Fasta(args['fasta2'])
    # choose program
    if args['program'] == 1:
        # export to fasta
        sys.stdout = open(args['output'], 'a')
        # iterate the following 2 lists
        for key in features1.keys():
            if key in features2.keys():
                print(''.join([">",features1[key].long_name]))
                print('\n'.join(split_every_60(features1[key][:].seq)))
        sys.stdout.close()
    else:
        # export to fasta
        sys.stdout = open(args['output'], 'a')
        # iterate the following 2 lists
        for key in features1.keys():
            if key not in features2.keys():
                print(''.join([">",features1[key].long_name]))
                print('\n'.join(split_every_60(features1[key][:].seq)))
        sys.stdout.close()

if __name__ == '__main__':
    main()