# python3
from gooey import *
from Bio import SeqIO
import os
# imput parameters
@Gooey(required_cols=5, program_name='genomic genbank to gene with annotations', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
	ap = GooeyParser(description="retrieve gene from genbank in genbank format containing the sequence, UTR's and exon-intron boundaries")
	ap.add_argument("-in", "--input", required=True, widget='FileChooser', help="input genomic genbank file")
	ap.add_argument("-chr", "--chr", required=True, type=str, help="chromosome/scaffold/contig the gene is located")
	ap.add_argument("-start", "--start", required=True, type=int, help="start of the gene in the chromosome/scaffold/contig")
	ap.add_argument("-end", "--end", required=True, type=int, help="end of the gene in the chromosome/scaffold/contig")
	ap.add_argument("-gb", "--genbank", required=True, type=str, widget='FileSaver', help="output genbank file")
	args = vars(ap.parse_args())
	# retrieve the gene with annotations from the genomic genbank file
	for record in SeqIO.parse(args['input'], "genbank"):
		if record.id == args['chr']:
			trimmed = record[int(args['start'] -1):args['end']]
	# export to genbank format
	SeqIO.write(trimmed, args['genbank'], "genbank")

if __name__ == '__main__':
	main()