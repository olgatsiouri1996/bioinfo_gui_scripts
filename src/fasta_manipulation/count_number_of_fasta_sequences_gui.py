# python 3
import os
from gooey import *
from pyfaidx import Fasta
# input parameters
@Gooey(required_cols=2, program_name='count number of fasta sequences', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser(description="count number of fasta sequences for each fasta file by indexing the  multi-fasta files in a user specified directory")
    ap.add_argument("-dir", "--directory", required=True, widget='DirChooser', help="directory of input mult-fasta files with extensions: .fasta, .fa, .fna, .faa, .fsa, .ffn, frn, .mpfa")
    ap.add_argument("-out", "--output", required=True, widget='FileSaver', help="output tab seperated 2-column txt file with fasta filenames and number of fasta sequences")
    args = vars(ap.parse_args())
    # create fasta index for each fasta file
    with open(args['output'], 'a') as filehandle:
        for filename in sorted(os.listdir(os.chdir(args['directory']))):
            if filename.endswith(".fa") or filename.endswith(".fasta") or filename.endswith(".fna") or filename.endswith(".faa") or filename.endswith(".fsa") or filename.endswith(".ffn") or filename.endswith(".frn") or filename.endswith(".mpfa"):
                features = Fasta(filename)
                filehandle.write('%s\n' % '\t'.join([filename,str(len(features.keys()))]))
                del features

if __name__ == '__main__':
    main()
