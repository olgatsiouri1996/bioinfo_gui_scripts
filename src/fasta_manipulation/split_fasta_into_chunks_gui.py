# python3
from gooey import *
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
# input parameters
@Gooey(required_cols=3, program_name= 'split fasta into chunks', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
	ap = GooeyParser()
	ap.add_argument("-in", "--input", required=True, widget='FileChooser', help="input single or multi-fasta file")
	ap.add_argument("-out", "--output", required=False, widget='FileSaver', help="output multi-fasta file")
	ap.add_argument("-step", "--step", required=True, help="step size for chunk creation, type = integer")
	ap.add_argument("-win", "--window", required=True, help="window size for chunk creation, type = integer")
	ap.add_argument("-pro", "--program", type=int, default=1, required=False, help="program output to select 1) 1 multi-fasta file 2) many single-fasta files")
	ap.add_argument("-dir", "--directory", required=False, type=str, widget='DirChooser', help="output directory to save the single-fasta files")
	args = vars(ap.parse_args())
	# main
	seqs = []
	headers = [] 
	chunks = []# setup empty lists
	# import multi or single-fasta file
	for record in SeqIO.parse(args['input'], "fasta"):
		for i in range(0, len(record.seq) - int(args['window']) + 1, int(args['step'])):
			seqs.append(record.seq[i:i + int(args['window'])])
			headers.append('_'.join([record.id,str(i + 1)]))
	# export to multi or single-fasta
	for (seq, header) in zip(seqs,headers):
		chunks.append(SeqRecord(Seq(str(seq)),id=str(header),description=""))
	if args['program']==1:
		SeqIO.write(chunks,args['output'], "fasta")
	else:
		# set working directory
		os.chdir(args['directory'])
		for chunk in chunks:	
			SeqIO.write(chunk, ''.join([str(chunk.id),".fasta"]), "fasta")

if __name__ == '__main__':
	main()
