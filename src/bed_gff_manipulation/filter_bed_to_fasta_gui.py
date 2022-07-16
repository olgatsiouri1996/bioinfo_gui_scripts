# python3
from gooey import *
import sys
from pyfaidx import Fasta
import pandas as pd
import warnings
# input parameters
@Gooey(required_cols=3, program_name= 'filter bed to fasta', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser()
    ap.add_argument("-bed", "--bed", required=True, widget='FileChooser', help="input bed file(made with bedops, every feature in the \'.gff\' or \'.gff3\' file should have an \'ID\' tag in the \'attributes\' column)")
    ap.add_argument("-in", "--input", required=True, widget='FileChooser', help="input fasta file")
    ap.add_argument("-out", "--output", required=True, widget='FileSaver', help="output fasta file")
    ap.add_argument("-fea", "--feature", required=False, default='gene',  type=str, help="specify the pattern to select the lines that have it")
    args = vars(ap.parse_args())
    # main
    # create function to split the input sequence based on a specific number of characters(60)
    def split_every_60(s): return [str(s)[i:i+60] for i in range(0,len(str(s)),60)]
    # ignore warnings
    warnings.filterwarnings('ignore')
    # import bed with no headers specified
    df = pd.read_csv(args['bed'], sep= "\t", header=None)
    # select the rows containing the feature
    bool2 = df.iloc[:, 7].str.contains(args['feature']) 
    df = df[bool2]
    # convert each column to list
    chrom = df.iloc[:,0].values.tolist()
    start = df.iloc[:,1].values.tolist()
    end = df.iloc[:,2].values.tolist()
    ids = df.iloc[:,3].values.tolist()
    strand = df.iloc[:,5].values.tolist()
    # import fasta file
    features = Fasta(args['input'])
    # iterate all below lists in pairs
    sys.stdout = open(args['output'], 'a')
    for (a, b, c, d, e) in zip(ids, chrom, start, end, strand):
        if str(e) == "+":
            print(''.join([">",str(a)," ",str(b),":",str(int(c) + 1),"-",str(d)]))
            print('\n'.join(split_every_60(features[str(b)][int(c):int(d)].seq)))
        else:
            print(''.join([">",str(a)," ",str(b),":",str(int(c) + 1),"-",str(d)," ","reverse complement"]))
            print('\n'.join(split_every_60(features[str(b)][int(c):int(d)].reverse.complement.seq)))
    sys.stdout.close()

if __name__ == '__main__':
    main()
