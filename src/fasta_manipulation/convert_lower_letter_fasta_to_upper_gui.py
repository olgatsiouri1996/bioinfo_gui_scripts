# python3
from gooey import *
from pyfaidx import Fasta
# input parameters
@Gooey(required_cols=1, program_name= 'convert lower letter fasta to upper', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser(description="ovewrite and convert a fasta with lowercase letters to uppercase")
    ap.add_argument("-in", "--input", required=True, widget='FileChooser', help="input fasta file")
    args = vars(ap.parse_args())
    # main
    # import fasta file
    features = Fasta(args['input'],mutable=True)
    # iterate all below lists in pairs
    for key in features.keys():
        features[key][:] = str(features[key][:].seq).upper()

if __name__ == '__main__':
    main()

