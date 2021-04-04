# python3
from gooey import *
from Bio import SeqIO
import sys
# imput parameters
@Gooey(required_cols=5, program_name='codon optimization', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser()
    ap.add_argument("-up", "--upstream", required=True, widget='FileChooser', help="upstream adapter(fasta format)")
    ap.add_argument("-fa", "--fasta", required=True, widget='FileChooser', help="input fasta file with the sequence to add adapters")
    ap.add_argument("-down", "--downstream", required=True, widget='FileChooser', help="downstream adapter(fasta format)")
    ap.add_argument("-id", "--seqid", required=True, help="fasta header of the output file")
    ap.add_argument("-out", "--output", required=True, widget='FileSaver', help="output fasta file")
    args = vars(ap.parse_args())
# main
# first adapter
    for record in SeqIO.parse(args['upstream'], "fasta"):
        x = str(record.seq)
# sequence
    for record in SeqIO.parse(args['fasta'], "fasta"):
        y = str(record.seq)
# second adapter
    for record in SeqIO.parse(args['downstream'], "fasta"):
        z = str(record.seq)
# merge
    seqad = x + y + z
# output
    sys.stdout = open(args['output'], 'a')
    print(">"+args['seqid'], seqad, sep='\n')
    sys.stdout.close()

if __name__ == '__main__':
    main()