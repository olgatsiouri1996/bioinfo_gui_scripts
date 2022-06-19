# python3
from gooey import *
from pyfaidx import Fasta
# input parameters
@Gooey(required_cols=2, program_name= 'fasta to bed', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser()
    ap.add_argument("-in", "--input", required=True, widget='FileChooser', help="input multi-fasta file")
    ap.add_argument("-bed", "--bed", required=True, widget='FileSaver', help="output bed file with id start and end locations")
    args = vars(ap.parse_args())
    # main
    # import fasta file
    features = Fasta(args['input'])
    # export to a bed file
    with open(args['bed'], 'w') as filehandle:
        for key in features.keys():
            filehandle.write('%s\n' % '\t'.join([key,"1",str(features[key][:].end)]))

if __name__ == '__main__':
    main()