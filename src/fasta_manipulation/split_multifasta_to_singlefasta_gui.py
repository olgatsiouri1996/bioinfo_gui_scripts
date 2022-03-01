# python3
from gooey import *
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
# imput parameters
@Gooey(required_cols=1, program_name= 'split multi-fasta to single-fasta', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
	ap = GooeyParser()
	ap.add_argument("-mfa", "--multifasta", required=True, widget='FileChooser', help="input multi-fasta file to split to single-fasta")
	ap.add_argument("-dir", "--directory", required=False, type=str, widget='DirChooser', help="output directory to save the single-fasta files")
	args = vars(ap.parse_args())
# main
	records = SeqIO.parse(args['multifasta'], "fasta")
# set working directory
	os.chdir(args['directory'])
	for record in records:
		SeqIO.write(record, ''.join([record.id,".fasta"]), "fasta")

if __name__ == '__main__':
	main()
