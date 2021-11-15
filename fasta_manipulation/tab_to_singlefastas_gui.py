# python3
import itertools
from gooey import *
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
# input arguments
@Gooey(required_cols=1, program_name='tabular to single-fasta files', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser(description="convert each row of a tabular file with the fasta headers and sequences in each row in single-fasta files")
    ap.add_argument("-in", "--input", required=True, widget="FileChooser", help="input txt file")
    args = vars(ap.parse_args())
	# main
    df = pd.read_csv(args['input'], header=None, sep="\t")
	# select ids and sequence columns, convert to lists
    headers = df.iloc[:,0].values.tolist()
    sequences = df.iloc[:,1].values.tolist()
	# iter elements on pairs to export to single-fasta files
    for (ids, seq) in zip(headers, sequences):
		    seq_for_fasta=SeqRecord(Seq(str(seq)),id=str(ids),description="")
		    SeqIO.write(seq_for_fasta, "".join([str(ids),".fasta"]), "fasta")

if __name__ == '__main__':
    main()