# python3
from gooey import *
from Bio import SeqIO
# imput parameters
@Gooey(required_cols=1, program_name=' fasta subset by length', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
	ap = GooeyParser(description="filters a fasta file by keeping sequences or headers under a range of length")
	ap.add_argument("-in", "--input", required=True,widget='FileChooser', help="input fasta file")
	ap.add_argument("-pro", "--program",type=int, default=1, required=False, help="choose to output a fasta file or a txt file with headers(1.fasta file with sequence length in fasta description, 2.txt file with headers).")
	ap.add_argument("-out", "--output", required=False, widget='FileSaver', help="output fasta file")
	ap.add_argument("-max", "--max", required=False, default=300, help="max number of sequence length")
	ap.add_argument("-min", "--min", required=False, default=1, help="min number of sequence length")
	ap.add_argument("-headers", "--headers", required=False, widget='FileSaver', help="file to save the output fasta headers")
	args = vars(ap.parse_args())
	# main
	# choose program
	program = args['program']
	# select sequences
	if program == 1:
	    sequences = []  # setup an empty list
	    for record in SeqIO.parse(args['input'], "fasta"):
	        if int(args['min']) <= len(record.seq) <= int(args['max']):
	            # add this record to the list
	            record.description = ''.join(["length:"," ",str(len(record.seq))," ",str(record.description).split(record.id)[1]])
	            sequences.append(record)
	    # export to fasta1
	    SeqIO.write(sequences, args['output'], "fasta")
	# retrieve headers only
	else:
	    headers = []  # setup an empty list
	    for record in SeqIO.parse(args['input'], "fasta"):
	        if int(args['min']) <= len(record.seq) <= int(args['max']):
	            # add this record to the list
	            headers.append(record.id)
	    # export to txt
	    with open(args['headers'], 'w') as filehandle:
	        for listitem in headers:
	            filehandle.write('%s\n' % listitem)

if __name__ == '__main__':
	main()