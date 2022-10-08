# python3
import os
from gooey import *
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
# input arguments
@Gooey(required_cols=2, program_name='tabular to single-fasta files', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser(description="convert each row of a tabular file with the fasta headers fasta descriptions and sequences in each row in single-fasta files")
    ap.add_argument("-in", "--input", required=True, widget="FileChooser", help="input txt file")
    ap.add_argument("-dir", "--directory", required=True, type=str, widget="DirChooser", help="directory of output fasta files")
    args = vars(ap.parse_args())
	# main
    df = pd.read_csv(args['input'], header=None, sep="\t")
    # select ids and sequence columns, convert to lists
    headers = df.iloc[:,0].values.tolist()
    descriptions = df.iloc[:,1].fillna("").values.tolist()   
    sequences = df.iloc[:,2].values.tolist()	
    # select output directory for single-fasta files
    os.chdir(args['directory'])
	# iter elements on pairs to export to single-fasta files
    for (ids, seq, desc) in zip(headers, sequences, descriptions):
        seq_for_fasta=SeqRecord(Seq(str(seq)),id=str(ids),description=str(desc))
        SeqIO.write(seq_for_fasta, "".join([str(ids),".fasta"]), "fasta")

if __name__ == '__main__':
    main()
