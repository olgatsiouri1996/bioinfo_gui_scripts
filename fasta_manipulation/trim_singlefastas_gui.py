# python3
import os
from gooey import *
from Bio import SeqIO
# input parameters
@Gooey(required_cols=3, program_name='trim multiple single-fasta files', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser()
    ap.add_argument("-start", "--start_fasta", required=True, type=int, help="region to start writing the fasta file(min number 0)")
    ap.add_argument("-stop", "--stop_fasta", required=True, type=int, help="region to stop writing the fasta file(negative number to remove nucleotides from the end of the sequence")
    ap.add_argument("-dir", "--directory", required=True, type=str, widget='DirChooser', help="directory to search for fasta files")
    args = vars(ap.parse_args())
    # main
    # import each fasta file from a working directory of choice
    for filename in sorted(os.listdir(os.chdir(args['directory']))):
        if filename.endswith(".fa") or filename.endswith(".fasta"):
            # read each file, trim and create SeqRecord to export
            record = SeqIO.read(filename, "fasta")
            sequence = record[args['start_fasta']:args['stop_fasta']]
            # export to fasta
            SeqIO.write(sequence, "".join([filename.split(".")[0],"_","trimmed",".fasta"]), "fasta")


if __name__ == '__main__':
    main()