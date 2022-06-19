# python3
from gooey import *
from pyfaidx import Fasta
# input parameters
@Gooey(required_cols=2, program_name= 'retrieve full fasta headers', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser(description="retrieve fasta identifiers and their descriptions")
    ap.add_argument("-in", "--input", required=True, widget='FileChooser', help="input multi-fasta file")
    ap.add_argument("-headers", "--headers", required=True, widget='FileSaver', help="1-column txt file to save the output fasta headers with descriptions")
    args = vars(ap.parse_args())
    # main
    # import fasta file
    features = Fasta(args['input'])
    # export to 1-column txt file
    with open(args['headers'], 'w') as filehandle:
        for key in features.keys():
            filehandle.write('%s\n' % features[str(key)].long_name)

if __name__ == '__main__':
    main()