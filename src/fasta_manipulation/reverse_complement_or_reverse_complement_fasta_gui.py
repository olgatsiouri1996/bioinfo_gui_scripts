# python3
from gooey import *
from Bio import SeqIO
# input parameters
@Gooey(required_cols=2, program_name='reverse, complement or reverse complement fasta', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser()
    ap.add_argument("-in", "--input", required=True, widget='FileChooser', help="input fasta file")
    ap.add_argument("-pro", "--program", required=False, default=1, type=int, help="program to choose 1. reverse complement, 2. reverse, 3. complement. Default is 1")
    ap.add_argument("-out", "--output", required=True, widget='FileSaver', help="output fasta file")
    args = vars(ap.parse_args())
    # main
    sequences = []  # setup an empty list
    # select program
    if args['program'] == 1:
        for record in SeqIO.parse(args['input'], "fasta"):
            # add this record to the list
            record.seq = record.seq.reverse_complement()
            sequences.append(record)
            SeqIO.write(sequences, args['output'], "fasta")
    elif args['program'] == 2:
        for record in SeqIO.parse(args['input'], "fasta"):
            # add this record to the list
            sequences.append(record[::-1])
            SeqIO.write(sequences, args['output'], "fasta")
    else:
        for record in SeqIO.parse(args['input'], "fasta"):
            # add this record to the list
            record.seq = record.seq.complement()
            sequences.append(record)
            SeqIO.write(sequences, args['output'], "fasta")

if __name__ == '__main__':
    main()