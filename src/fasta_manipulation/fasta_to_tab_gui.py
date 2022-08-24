# python3
from gooey import *
import sys
import pyfastx
# input parameters
@Gooey(required_cols=2, program_name='fasta to tab', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser()
    ap.add_argument("-in","--input", required=True, widget='FileChooser',help="input multi-fasta file(all fasta records should either have no fasta description or all of them should have fasta description)")
    ap.add_argument("-pro", "--program", required=False, default=1, type=int, help="output to choose: 1) 2-column txt file with fasta identifiers and fasta sequences 2) 2-column txt file with full fasta headers and fasta sequences 3) 3-column txt file with fasta identifiers, fasta descriptions and fasta sequences")
    ap.add_argument("-out","--output", required=True, widget='FileSaver', help="output txt file")
    args = vars(ap.parse_args())
    # main
    # choose program
    if args['program'] == 1:
        sys.stdout = open(args['output'], 'w')
        for name,seq,comment in pyfastx.Fastx(args['input']):
            print(name,seq,sep='\t')  
        sys.stdout.close()
    elif args['program'] == 2:
        sys.stdout = open(args['output'], 'w')
        for name,seq,comment in pyfastx.Fastx(args['input']):
            print(' '.join([name,comment]),seq,sep='\t')
        sys.stdout.close()
    else:
        sys.stdout = open(args['output'], 'w')
        for name,seq,comment in pyfastx.Fastx(args['input']):
            print(name,comment,seq,sep='\t')  
        sys.stdout.close()

if __name__ == '__main__':
    main()