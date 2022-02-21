# python3
from gooey import *
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
# imput parameters
@Gooey(required_cols=1, program_name= 'split multi-fasta to single-fasta', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
	ap = GooeyParser()
	ap.add_argument("-mfa", "--multifasta", required=True, widget='FileChooser', help="input multi-fasta file to split to single-fasta")
	args = vars(ap.parse_args())
	# main
	for record in SeqIO.parse(args['multifasta'], "fasta"):
		one_seq = SeqRecord(Seq(record.seq),id=record.id,description="")
		SeqIO.write(one_seq, ''.join([record.id,".fasta"]), "fasta")

if __name__ == '__main__':
	main()
