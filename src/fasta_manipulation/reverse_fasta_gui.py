# python3
from gooey import *
from Bio import SeqIO
# input parameters
@Gooey(required_cols=2, program_name='reverse fasta', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser(description="reverse the sequences of a multi or single-fasta file")
    ap.add_argument("-in", "--input", required=True, widget='FileChooser', help="input fasta file")
    ap.add_argument("-out", "--output", required=True, widget='FileSaver', help="output fasta file")
    args = vars(ap.parse_args())
    # main
    sequences = []  # setup an empty list
    for record in SeqIO.parse(args['input'], "fasta"):
        # add this record to the list
        sequences.append(record[::-1])
    SeqIO.write(sequences, args['output'], "fasta")

if __name__ == '__main__':
    main()