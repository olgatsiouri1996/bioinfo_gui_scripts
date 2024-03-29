# python3
from gooey import *
import os
import sys
from pyfaidx import Fasta
# imput parameters
@Gooey(required_cols=2, program_name= 'split multi-fasta to single-fasta', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
	ap = GooeyParser()
	ap.add_argument("-mfa", "--multifasta", required=True, widget='FileChooser', help="input multi-fasta file to split to single-fasta")
	ap.add_argument("-dir", "--directory", required=True, type=str, widget='DirChooser', help="output directory to save the single-fasta files")
	args = vars(ap.parse_args())
# create function to split the input sequence based on a specific number of characters(60)
	def split_every_60(s): return [str(s)[i:i+60] for i in range(0,len(str(s)),60)]
# import multi-fasta file
	features = Fasta(args['multifasta'])
# select output directory
	os.chdir(args['directory'])
# export each record to a single-fasta    
	for key in features.keys():
	    sys.stdout = open(''.join([str(key),".fasta"]), 'a')
	    print(''.join([">",features[str(key)].long_name]).replace('\r',''))
	    print('\n'.join(split_every_60(features[str(key)][:].seq)))
	    sys.stdout.close()

if __name__ == '__main__':
	main()
