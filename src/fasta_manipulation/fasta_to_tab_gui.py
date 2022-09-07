# python3
from gooey import *
import sys
from pyfaidx import Fasta
# input parameters
@Gooey(required_cols=2, program_name='fasta to tab', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser()
    ap.add_argument("-in","--input", required=True, widget='FileChooser',help="input multi-fasta file")
    ap.add_argument("-pro", "--program", required=False, default=1, type=int, help="output to choose: 1) 2-column txt file with fasta identifiers and fasta sequences 2) 2-column txt file with full fasta headers and fasta sequences 3) 3-column txt file with fasta identifiers, fasta descriptions and fasta sequences")
    ap.add_argument("-out","--output", required=True, widget='FileSaver', help="output txt file")
    args = vars(ap.parse_args())
    # main
    # index multi-fasta file
    features = Fasta(args['input'])
    # choose program
    if args['program'] == 1:
        sys.stdout = open(args['output'], 'w')
        for key in features.keys():
            print(key,features[key][:].seq,sep='\t')  
        sys.stdout.close()
    elif args['program'] == 2:
        sys.stdout = open(args['output'], 'w')
        for key in features.keys():
            print(features[key].long_name,features[key][:].seq,sep='\t')
        sys.stdout.close()
    else:
        sys.stdout = open(args['output'], 'w')
        for key in features.keys():
            try:
                description = str(features[key].long_name).split(' ',1)[1]
            except IndexError:
                description = ""
     
            print(key,description,features[key][:].seq,sep='\t')  
        sys.stdout.close()

if __name__ == '__main__':
    main()
