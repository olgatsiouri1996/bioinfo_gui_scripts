# python3
from gooey import * 
import os
from Bio import SeqIO
import pandas as pd
# imput parameters
@Gooey(required_cols=2, program_name='single-fastas to tabular txt file', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
	ap = GooeyParser()
	ap.add_argument("-dir", "--directory", required=True, type=str, widget='DirChooser', help="directory to search for fasta files")
	ap.add_argument("-out","--output", required=True, widget='FileSaver', help="output txt file with 3 columns without column names consisting of: id, description, sequence")
	args = vars(ap.parse_args())
# main
	seqs = []
	ids = [] 
	descriptions = [] # setup empty lists
# import each fasta file from the working directory
	for filename in sorted(os.listdir(os.chdir(args['directory']))):
		if filename.endswith(".fa") or filename.endswith(".fasta"):
			for record in SeqIO.parse(filename, "fasta"):
				ids.append(record.id)
				try:
					descriptions.append(record.description.split(' ',1)[1])
				except IndexError:
					descriptions.append('')
				seqs.append(record.seq)
# put the 2 list in a data frame of 2 columns
	dfasta = pd.DataFrame()
	dfasta['id'] = ids
	dfasta['desc'] = descriptions
	dfasta['seq'] = seqs
# export data frame to a tabular txt file
	with open(args['output'], 'a') as f:
	    f.write(
	        dfasta.to_csv(header = False, index = False, sep= "\t", lineterminator= '\n')
	    )

if __name__ == '__main__':
    main()
