# python3
from gooey import *
from pyfaidx import Fasta
# input arguments
@Gooey(required_cols=2, program_name='fasta to header', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser(description="exctract the headers from a multifasta file")
    ap.add_argument("-in", "--input", required=True, widget='FileChooser', help="input fasta file")
    ap.add_argument("-headers", "--headers", required=True, widget='FileSaver', help="1-column txt file to save the output fasta headers")
    args = vars(ap.parse_args())
# main
# index multi-fasta file
    features = Fasta(args['input'])
# export to 1-column txt file
    with open(args['headers'], 'w') as filehandle:
        for key in features.keys():
            filehandle.write('%s\n' % key)

if __name__ == '__main__':
    main()
