# python3
from gooey import *
from pyfaidx import Fasta
# input parameters
@Gooey(required_cols=2, program_name= 'softmask fasta with bed', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser(description="ovewrite and softmask a multi-fasta file")
    ap.add_argument("-bed", "--bed", required=True, widget='FileChooser', help="input bed file(made with bedops)")
    ap.add_argument("-in", "--input", required=True, widget='FileChooser', help="input fasta file")
    args = vars(ap.parse_args())
    # main
    # setup empty lists
    chrom = []
    start = []
    end = []
    # import bed 
    with open(args['bed'], 'r') as f:
        for line in f:
            # convert each column to list
            chrom.append(line.split()[0])
            start.append(line.split()[1])
            end.append(line.split()[2])
    # import fasta file
    features = Fasta(args['input'],mutable=True)
    # iterate all below lists in pairs
    for (a, b, c) in zip(chrom, start, end):
        features[str(a)][int(b):int(c)] = str(features[str(a)][int(b):int(c)].seq).lower()

if __name__ == '__main__':
    main()

