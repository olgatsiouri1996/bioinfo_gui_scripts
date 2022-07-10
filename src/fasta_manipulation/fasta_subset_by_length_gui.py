# python3
from gooey import *
from pyfaidx import Fasta
import sys
# imput parameters
@Gooey(required_cols=1, program_name=' fasta subset by length', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser(description="filters a fasta file by keeping sequences or headers under a range of length")
    ap.add_argument("-in", "--input", required=True,widget='FileChooser', help="input fasta file")
    ap.add_argument("-pro", "--program",type=int, default=1, required=False, help="choose to output a fasta file or a txt file with headers(1.fasta file with sequence length in fasta description, 2.txt file with headers).")
    ap.add_argument("-out", "--output", required=False, widget='FileSaver', help="output fasta file")
    ap.add_argument("-max", "--max", required=False, default=300, help="max number of sequence length")
    ap.add_argument("-min", "--min", required=False, default=1, help="min number of sequence length")
    ap.add_argument("-headers", "--headers", required=False, widget='FileSaver', help="file to save the output fasta headers")
    args = vars(ap.parse_args())
    # main
    # create function to split the input sequence based on a specific number of characters(60)
    def split_every_60(s): return [str(s)[i:i+60] for i in range(0,len(str(s)),60)]
    # create fasta index
    features = Fasta(args['input'])
    # select sequences
    # choose program
    if args['program'] == 1:
        sys.stdout = open(args['output'], 'a')
        for key in features.keys():
            if int(args['min']) <= features[key][:].end <= int(args['max']):
                print(''.join([">",features[key].long_name," ","length:"," ",str(features[key][:].end)]))
                print('\n'.join(split_every_60(features[key][:].seq)))
        sys.stdout.close()
    # retrieve headers only
    else:
        # export to txt
        with open(args['headers'], 'w') as filehandle:
            for key in features.keys():
                if int(args['min']) <= features[key][:].end <= int(args['max']):
                    filehandle.write('%s\n' % key)

if __name__ == '__main__':
    	main()
