# python 3
import os
import glob
from gooey import *
from pyfaidx import Fasta
# input parameters
@Gooey(required_cols=2, program_name='count number of fasta sequences', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser(description="count number of fasta sequences for each fasta file by indexing the  multi-fasta files in a user specified directory")
    ap.add_argument("-dir", "--directory", required=True, widget='DirChooser', help="directory of input mult-fasta files with extensions: .fasta, .fa, .fna, .faa, .fsa, .ffn, frn, .mpfa")
    ap.add_argument("-out", "--output", required=True, widget='FileSaver', help="output tab seperated 2-column txt file with fasta filenames and number of fasta sequences")
    ap.add_argument("-type", "--type", required=False, type=str, default="all fasta extensions", choices=["all fasta extensions", "a specific pattern","many patterns"],  help="input pattern type")
    ap.add_argument("-pat", "--pattern", required=False, type=str,  help="input pattern of files to count the sequences for")
    ap.add_argument("-pats", "--patterns", required=False, widget='FileChooser',  help="1-column txt file with patterns")
    args = vars(ap.parse_args())
    # choose programs
    if args['type']=="all fasta extensions":
        filepaths = tuple(glob.glob(os.path.join(args['directory'], "*.fa")) + glob.glob(os.path.join(args['directory'], "*.fasta")) + glob.glob(os.path.join(args['directory'], "*.fna")) + glob.glob(os.path.join(args['directory'], "*.ffn")) + glob.glob(os.path.join(args['directory'], "*.faa")) + glob.glob(os.path.join(args['directory'], "*.frn")))
    elif args['type']=="a specific pattern":
        filepaths = tuple(glob.glob(os.path.join(args['directory'], args['pattern'])))
    else:
        filepaths = set()
        file_patters = (str(line.rstrip()) for line in open(args['patterns']))
        for file_pattern in file_patters:
            filepaths.update(glob.glob(os.path.join(args['directory'], file_pattern)))
    
    def count_seqs(filepath):
        filename = os.path.splitext(os.path.basename(filepath))[0]
        fasta_file = Fasta(filepath)
        return f"{filename}\t{str(len(fasta_file.keys()))}"

    counted_data = map(count_seqs,filepaths)

    with open(args['output'],'w') as output_txt:
        output_txt.write(f'{"filename"}\t{"number_of_sequences"}\n')
        output_txt.write("\n".join(counted_data))

        # Remove .fai files in the input directory
        for fai_file in glob.glob(os.path.join(args['directory'], "*.fai")):
            os.remove(fai_file)

if __name__ == '__main__':
    main()
