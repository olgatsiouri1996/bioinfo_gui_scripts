# python3
from gooey import *
from Bio import SeqIO
import pandas as pd
import os
# input parameters
@Gooey(required_cols=4, program_name='kmer generator windows', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser(description="splits a fasta file with user specified length and fragment overlap")
    ap.add_argument("-in", "--input", required=True, widget='FileChooser', help="input fasta file")
    ap.add_argument("-step", "--step", required=True, help="step size to split fasta, type = int")
    ap.add_argument("-win", "--window", required=True, help="window size of splitted subsets, type = int")
    ap.add_argument("-out", "--output", required=True, widget='FileSaver', help="output fasta file")
    args = vars(ap.parse_args())
# main
    sequences = []
    headers = [] # setup empty lists
    for record in SeqIO.parse(args['input'], "fasta"):
        for i in range(0, len(record.seq) - int(args['window']) + 1, int(args['step'])):
            sequences.append(record.seq[i:i + int(args['window'])])
            headers.append(i)
# create data frame
    df = pd.DataFrame()
    df['id'] = headers
    df['seq'] = sequences
# export
    with open("out.tab", 'a') as f:
        f.write(
            df.to_csv(header = False, index = False, sep = '\t', doublequote= False, line_terminator= '\n')
        )

# convert to fasta
    convert = SeqIO.convert("out.tab", "tab", args['output'], "fasta")
    os.system("del out.tab")

if __name__ == '__main__':
    main()
