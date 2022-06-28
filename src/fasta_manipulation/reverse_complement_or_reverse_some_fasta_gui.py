# python3
from gooey import *
from Bio import SeqIO
# input parameters
@Gooey(required_cols=3, program_name='reverse complement or reverse some fasta', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser(description="reverse complement or reverse some sequences in a multi-fasta file")
    ap.add_argument("-in", "--input", required=True, widget='FileChooser', help="input fasta file")
    ap.add_argument("-ids", "--ids", required=True, widget='FileChooser', help="file with fasta headers to reorient some output fasta sequences")
    ap.add_argument("-pro", "--program", required=False, default=1, type=int, help="program to choose 1. reverse complement, 2. reverse")
    ap.add_argument("-out", "--output", required=True, widget='FileSaver', help="output fasta file")
    args = vars(ap.parse_args())
# main
# import the txt file with headers you want to reorient the sequence from the input multi-fasta
    with open(args['ids'], 'r') as f:
        headers = f.readlines()
    headers = [x.strip() for x in headers]
# setup an empty list
    sequences = [] 
# choose program
    if args['program'] == 1:
        for record in SeqIO.parse(args['input'], "fasta"):
            if record.id in headers:
            # add this record to the list
                record.seq = record.seq.reverse_complement()
                sequences.append(record)
            else:
                sequences.append(record)
# export to fasta
        SeqIO.write(sequences, args['output'], "fasta")
    else:
        for record in SeqIO.parse(args['input'], "fasta"):
            if record.id in headers:
            # add this record to the list
                sequences.append(record[::-1])
            else:
                sequences.append(record)
# export to fasta
        SeqIO.write(sequences, args['output'], "fasta")

if __name__ == '__main__':
    main()
