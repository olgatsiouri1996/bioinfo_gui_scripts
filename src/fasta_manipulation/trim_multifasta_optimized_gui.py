# python3
from gooey import *
import sys
from pyfaidx import Fasta
# input parameters
@Gooey(required_cols=2, program_name= 'trim multi-fasta optimized', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser(description="index and trim an imput multi-fasta file")
    ap.add_argument("-in", "--input", required=True, widget='FileChooser', help="input multi-fasta file")
    ap.add_argument("-start", "--start", required=False, default=0, type=int, help="number of nucleotides/amino acids to trim from the start of the sequence")
    ap.add_argument("-end", "--end", required=False, default=0, type=int, help="number of nucleotides/amino acids to trim from the start of the sequence")
    ap.add_argument("-out", "--output", required=True, widget='FileSaver', help="output multi-fasta file")
    args = vars(ap.parse_args())
    # main
    # create function to split the input sequence based on a specific number of characters(60)
    def split_every_60(s): return [str(s)[i:i+60] for i in range(0,len(str(s)),60)]
    # import fasta file
    features = Fasta(args['input'])
    # export trimmed sequences to an output multi-fasta file
    sys.stdout = open(args['output'], 'a')
    for key in features.keys():
        print(''.join([">",features[key].long_name," ",str(features[str(key)][args['start']:int(features[key][:].end - args['end'])].start),"-",str(features[str(key)][args['start']:int(features[key][:].end - args['end'])].end)]))
        print('\n'.join(split_every_60(features[str(key)][args['start']:int(features[key][:].end - args['end'])].seq)))
    sys.stdout.close()

if __name__ == '__main__':
    main()
