# python3
import sys
from gooey import *
from pyfaidx import Fasta
import re
# input parameters
@Gooey(required_cols=2, program_name='sort multi-fasta ascendingly', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser(description="sorts a multi-fasta with fasta identifiers(which after the prefix do not have leading 0s), consisting of letters and integers ascendingly, e.g. foo1, foo200, foo5000")
    ap.add_argument("-in", "--input", required=True, widget='FileChooser', help="input multi-fasta file")
    ap.add_argument("-out", "--output", required=True, widget='FileSaver', help="output multi-fasta file")
    args = vars(ap.parse_args())
# main
# create function to split the input sequence based on a specific number of characters(60)
    def split_every_60(s): return [str(s)[i:i+60] for i in range(0,len(str(s)),60)]
# import fasta file
    features = Fasta(args['input'])
    headers = [] # setup empty list
# create pattern to split
    r = re.compile("([a-zA-Z]+)([0-9]+)")
# retrieve prefix
    prefix = r.match(features[0].name).group(1)
# split prefix from fasta ids
    for key in features.keys():
        headers.append(int(key.split(prefix)[1]))
# sort in ascending order
    headers.sort()
# add removed prefix
    sorted_headers = [prefix + str(x) for x in headers]
# iterate sorted headers to extract sequences and export as multi-fasta
    sys.stdout = open(args['output'], 'a')
    for header in sorted_headers:
        print(''.join([">",features[str(header)].long_name]))
        print('\n'.join(split_every_60(features[str(header)][:].seq)))
    sys.stdout.close()

if __name__ == '__main__':
    main()