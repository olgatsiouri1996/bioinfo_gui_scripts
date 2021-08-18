# python3
from gooey import *
from Bio import SeqIO
# input parameters
@Gooey(required_cols=4, program_name='select headers by %aa content', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
	ap = GooeyParser(description="select the headers or the sequences from a fasta file by aminoacid content")
	ap.add_argument("-in", "--input", required=True, widget='FileChooser', help="input fasta file")
	ap.add_argument("-pro", "--program",type=int, default=1, required=False, help="choose to output a fasta file or a txt file with headers(1.fasta, 2.txt file with headers. Defaults to 1)")
	ap.add_argument("-max", "--max", required=True, help="max threshold of aa content, type = float")
	ap.add_argument("-min", "--min", required=True, help="min threshold of aa content, type = float")
	ap.add_argument("-aa", "--aa", required=True, help="aa to search the content for")
	ap.add_argument("-headers", "--headers", required=False, widget='FileSaver', help="file to save the output fasta headers")
	ap.add_argument("-out", "--output", required=False, widget='FileSaver', help="output fasta file")
	args = vars(ap.parse_args())
	# create aa_content function
	def aa_content(seq):
		return round((seq.count(args['aa']) / len(seq)) * 100, 2)
	# main
	# choose program
	program = args['program']
	# select sequences
	if program == 1:
		sequences = []  # setup an empty list
		for record in SeqIO.parse(args['input'], "fasta"):
			if float(args['min']) < aa_content(record.seq) < float(args['max']):
				# add this record to the list
				sequences.append(record)
		# export to fasta
		SeqIO.write(sequences, args['output'], "fasta")
	# retrieve headers only
	elif program == 2:
		headers = []  # setup an empty list
		for record in SeqIO.parse(args['input'], "fasta"):
			if float(args['min']) < aa_content(record.seq) < float(args['max']):
				# add this record to the list
				headers.append(record.id)
		# export to txt
		with open(args['headers'], 'w') as filehandle:
			for listitem in headers:
				filehandle.write('%s\n' % listitem)

if __name__ == '__main__':
	main()
