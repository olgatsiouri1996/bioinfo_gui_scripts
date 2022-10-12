# python3
from gooey import *
from Bio import SeqIO
import os
import pandas as pd
# imput parameters
@Gooey(required_cols=1, program_name='genomic genbank to gene/s with annotations', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
	ap = GooeyParser(description="retrieve 1 gene or many from genbank, in genbank format containing the sequence, UTR's and exon-intron boundaries")
	ap.add_argument("-in", "--input", required=True, widget='FileChooser', help="input genomic genbank file")
	ap.add_argument("-txt", "--txt", required=False, widget='FileChooser', help="4-column txt file with chromosome/scaffold/contig names, start, end coordinates and gene ids/names")
	ap.add_argument("-num", "--number", required=False, type=str, default='one', widget='Dropdown', choices=['one','many'], help="number of genbank files to create")
	ap.add_argument("-chr", "--chr", required=False, type=str, help="chromosome/scaffold/contig the gene is located")
	ap.add_argument("-start", "--start", required=False, type=int, help="start of the gene in the chromosome/scaffold/contig")
	ap.add_argument("-end", "--end", required=False, type=int, help="end of the gene in the chromosome/scaffold/contig")
	ap.add_argument("-gb", "--genbank", required=False, widget='FileSaver', help="output genbank file")
	ap.add_argument("-dir", "--directory", required=False, type=str, widget='DirChooser', help="director to save output genbank files")
	args = vars(ap.parse_args())
	# choose number of output files
	if args['number']=='one':
		# retrieve the gene with annotations from the genomic genbank file
		for record in SeqIO.parse(args['input'], "genbank"):
			if record.id == args['chr']:
				trimmed = record[int(args['start'] -1):args['end']]
		# export to genbank format
		SeqIO.write(trimmed, args['genbank'], "genbank")
	else:
		trimmed = [] # setup empty list
		# import txt file
		df = pd.read_csv(args['txt'], header=None, sep="\t")
		chrom = df.iloc[:,0].values.tolist()
		start = df.iloc[:,1].values.tolist()
		end = df.iloc[:,2].values.tolist()
		identifier = df.iloc[:,3].values.tolist()
		# retrieve the gene with annotations from the genomic genbank file
		for (a, b, c) in zip(chrom, start, end):
			for record in SeqIO.parse(args['input'], "genbank"):
				if record.id == str(a):
					trimmed.append(record[int(int(b) -1):int(c)])
		# change output directory
		os.chdir(args['directory'])
		for (rec, ids) in zip(trimmed, identifier):
		# export to genbank format
			SeqIO.write(rec, ''.join([str(ids),'.gb']), "genbank")			

if __name__ == '__main__':
	main()
