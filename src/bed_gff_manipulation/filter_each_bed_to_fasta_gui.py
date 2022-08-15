# python3
from gooey import *
import sys
from pyfaidx import Fasta
import pandas as pd
import warnings
# input parameters
@Gooey(required_cols=2, program_name= 'filter each bed to fasta', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser()
    ap.add_argument("-bed", "--bed", required=True, widget='FileChooser', help="txt file with the base name of fasta with chromosomes or scaffolds and  bed files(made with bedops, every feature in the \'.gff\' or \'.gff3\' file should have an \'ID\' tag in the \'attributes\' column), to import from the current directory")
    ap.add_argument("-fea", "--feature", required=False, widget='FileSaver', default='gene',  type=str, help="specify the pattern to select the lines that have it. Default is \'gene\'")
    ap.add_argument("-dir", "--directory", required=True, type=str, widget='DirChooser',  help="directory with bed files")
    args = vars(ap.parse_args())
    # main
    # create function to split the input sequence based on a specific number of characters(60)
    def split_every_60(s): return [str(s)[i:i+60] for i in range(0,len(str(s)),60)]
    # ignore warnings
    warnings.filterwarnings('ignore')
    # import the txt file with filenames
    with open(args['bed'], 'r') as f:
        orgs = f.readlines()
    orgs = [x.strip() for x in orgs]
    # go to the directory the .bed files are located
    os.chdir(args['directory'])
    # iterate for each file in the current directory
    for org in orgs:
        # import bed with no headers specified
        df = pd.read_csv(''.join([org,".bed"]), sep= "\t", header=None)
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
        features = Fasta(''.join([org,".fasta"]))
        # iterate all below lists in pairs
        sys.stdout = open(''.join([org,"_",args['feature'],".fasta"]), 'a')
        for (a, b, c, d, e) in zip(ids, chrom, start, end, strand):
            if str(e) == "+":
                print(''.join([">",str(a)," ",str(b),":",str(int(c) + 1),"-",str(d)]).replace('\r', ''))
                print('\n'.join(split_every_60(features[str(b)][int(c):int(d)].seq)))
            else:
                print(''.join([">",str(a)," ",str(b),":",str(int(c) + 1),"-",str(d)," ","reverse complement"]).replace('\r', ''))
                print('\n'.join(split_every_60(features[str(b)][int(c):int(d)].reverse.complement.seq)))
        sys.stdout.close()

        del df; del bool2; del chrom; del start; del end; del ids; del strand; del features

if __name__ == '__main__':
    main()
