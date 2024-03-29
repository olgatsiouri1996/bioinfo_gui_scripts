# python3
from gooey import *
import sys
from pyfaidx import Fasta
# input parameters
@Gooey(required_cols=4, program_name= 'split fasta into chunks optimized', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser()
    ap.add_argument("-in", "--input", required=True, widget='FileChooser', help="input single or multi-fasta file")
    ap.add_argument("-out", "--output", required=True, widget='FileSaver', help="output multi-fasta file")
    ap.add_argument("-step", "--step", required=True, type=int, help="step size for chunk creation, type = integer")
    ap.add_argument("-win", "--window", required=True, type=int, help="window size for chunk creation, type = integer")
    args = vars(ap.parse_args())
    # main
    # create function to split the input sequence based on a specific number of characters(60)
    def split_every_60(s): return [str(s)[i:i+60] for i in range(0,len(str(s)),60)]
    # import multi or single-fasta file
    features = Fasta(args['input'])
    # iterate input headers to extract sequences and export as multi-fasta
    sys.stdout = open(args['output'], 'a')
    for key in features.keys():
        # slice based on window and step
        for i in range(0, features[key][:].end - args['window'] + 1, args['step']):
            print(''.join([">",key,"_start_",str(i+1),"_stop_",str(i + args['window'])]).replace('\r',''))
            print('\n'.join(split_every_60(features[key][i:i + args['window']].seq)))
    sys.stdout.close()

if __name__ == '__main__':
	main()
